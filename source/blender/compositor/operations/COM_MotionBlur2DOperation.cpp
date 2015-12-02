/*
 * Copyright 2011, Blender Foundation.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * Contributor: 
 *		Tod Baudais
 */

#include "COM_MotionBlur2DOperation.h"
#include "MEM_guardedalloc.h"
#include <iostream>
#include <limits>

#define INDEX_COL(x,y) ((y * getWidth() + x) * COM_NUM_CHANNELS_COLOR)
#define INDEX_VEC(x,y) ((y * getWidth() + x) * COM_NUM_CHANNELS_VECTOR)
#define INDEX_VAL(x,y) ((y * getWidth() + x) * COM_NUM_CHANNELS_VALUE)

#define MAX_SAMPLES 1024*16
#define DEBUG_RENDER_ALPHA_MASK 0

MotionBlur2DOperation::MotionBlur2DOperation() : NodeOperation()
{
	addInputSocket(COM_DT_COLOR);
	addInputSocket(COM_DT_COLOR);
	addInputSocket(COM_DT_VALUE); // ZBUF
	addInputSocket(COM_DT_VALUE); // OBJID
	addOutputSocket(COM_DT_COLOR);
	m_settings = NULL;
	m_inputImageProgram = NULL;
	m_inputSpeedProgram = NULL;
	m_inputDepthProgram = NULL;
	m_inputObjIDProgram = NULL;
	m_cachedInstance = NULL;
	setComplex(true);
}
void MotionBlur2DOperation::initExecution()
{
	initMutex();
	m_inputImageProgram = getInputSocketReader(0);
	m_inputSpeedProgram = getInputSocketReader(1);
	m_inputDepthProgram = getInputSocketReader(2);
	m_inputObjIDProgram = getInputSocketReader(3);
	m_cachedInstance = NULL;
}

void MotionBlur2DOperation::deinitExecution()
{
	deinitMutex();
	if (m_cachedInstance) {
		MEM_freeN(m_cachedInstance);
		m_cachedInstance = NULL;
	}
}

void *MotionBlur2DOperation::initializeTileData(rcti *rect)
{
	if (m_cachedInstance) {
		return m_cachedInstance;
	}
	
	lockMutex();
	if (m_cachedInstance == NULL) {
        MemoryBuffer *color = (MemoryBuffer *)m_inputImageProgram->initializeTileData(rect);
        MemoryBuffer *speed = (MemoryBuffer *)m_inputSpeedProgram->initializeTileData(rect);
        float *data = (float *)MEM_callocN(MEM_allocN_len(color->getBuffer()), "motion2D data buffer");

        if (m_settings->deep_mode) {
            MemoryBuffer *depth = (MemoryBuffer *)m_inputDepthProgram->initializeTileData(rect);
            MemoryBuffer *objid = (MemoryBuffer *)m_inputObjIDProgram->initializeTileData(rect);
            generateMotionBlurDeep(data, color, speed, depth, objid);
        } else {
            generateMotionBlur(data, color, speed);
        }

		m_cachedInstance = data;
	}
	unlockMutex();
	return m_cachedInstance;
}

void MotionBlur2DOperation::executePixel(float output[4], int x, int y, void *data)
{
	float *buffer = (float *) data;
	copy_v4_v4(output, &buffer[INDEX_COL(x,y)]);
}

bool MotionBlur2DOperation::determineDependingAreaOfInterest(rcti *input, ReadBufferOperation *readOperation, rcti *output)
{
	if (m_cachedInstance == NULL) {
		rcti newInput;
		newInput.xmax = getWidth();
		newInput.xmin = 0;
		newInput.ymax = getHeight();
		newInput.ymin = 0;
		return NodeOperation::determineDependingAreaOfInterest(&newInput, readOperation, output);
	}
	else {
		return false;
	}
}

// Bressenham line
void MotionBlur2DOperation::line (int x0, int y0, int x1, int y1, Sample *samples, int *num_samples) {

    int dx = abs(x1-x0), sx = x0<x1 ? 1 : -1;
    int dy = abs(y1-y0), sy = y0<y1 ? 1 : -1; 
    int err = (dx>dy ? dx : -dy)/2, e2;
    *num_samples = 0;

    for(;;){
        samples[*num_samples].x = x0;
        samples[*num_samples].y = y0;
        ++(*num_samples);
        if (*num_samples >= MAX_SAMPLES)    // Too many samples
            break;

        if (x0==x1 && y0==y1)
            break;
        e2 = err;
        if (e2 >-dx) {  err -= dy; x0 += sx;    }
        if (e2 < dy) {  err += dx; y0 += sy;    }
    }
}

void MotionBlur2DOperation::generateMotionBlur(float *data, MemoryBuffer *color, MemoryBuffer *speed)
{

    // Allocate a temporary mask
    float *mask = (float *)MEM_callocN(MEM_allocN_len(color->getBuffer()), "motion2D mask buffer");

    // Clamp Multisample settings
    int multisample = m_settings->multisample;
    if (multisample < 0)        multisample = 1;
    else if (multisample > 8)   multisample = 8;
    int multisample2 = multisample * multisample;

    // Adjusted width and height given multisampling
    int width_multisample = getWidth() * multisample;
    int height_multisample = getHeight() * multisample;

    Sample line_samples[MAX_SAMPLES];
    int num_line_samples;

    // Blurring of pixels
    //
    // All this does is blurs the pixel contribution along the direction of movement for the pixel. If
    // multisampling is used then it acts internally like the images are X pixels larger (so
    // this gets slower) while taking care to keep the reading and writing of buffers correctly.
    //

    float forward_factor = m_settings->amount * multisample;
    float backwards_factor = m_settings->amount * multisample;

    for (int y = 0; y < height_multisample; ++y) {
        for (int x = 0; x < width_multisample; ++x) {
            int xm = x/multisample;
            int ym = y/multisample;

            // Render sample
            int index_col = INDEX_COL(xm,ym);
            float *color_pixel = color->getBuffer() + index_col;

            // Skip full transparent pixels
            if (color_pixel[3] == 0.0f)
                continue;

            int index_vec = INDEX_COL(xm,ym);
            float *speed_pixel = speed->getBuffer() + index_vec;
            float *mask_pixel = mask + index_col;

            // Calculate a motion vector
            int x0 = x + speed_pixel[0] * forward_factor;
            int y0 = y + speed_pixel[1] * forward_factor;
            int x1 = x - speed_pixel[0] * backwards_factor;
            int y1 = y - speed_pixel[1] * backwards_factor;

            // Build blur samples
            line (x0, y0, x1, y1, line_samples, &num_line_samples);
            float weight = 1.0f / (num_line_samples * multisample2);

            // Alpha hole mask
            float alpha_mask = 0.0f;

            for (int s = 0; s < num_line_samples; ++s) {
                int xs = line_samples[s].x;
                int ys = line_samples[s].y;

                // If outside image, assume alpha = 1
                if (xs < 0 || xs >= width_multisample || ys < 0 || ys >= height_multisample) {
                    alpha_mask += weight;

                } else {
                    int index_s = INDEX_COL(xs/multisample,ys/multisample);

                    // Read source alpha
                    float *alpha_mask_src_pixel = color->getBuffer() + index_s;
                    alpha_mask += alpha_mask_src_pixel[3] * weight;

                    // Write destination pixel
                    float *data_pixel = data + index_s;
                    data_pixel[0] += color_pixel[0] * weight;
                    data_pixel[1] += color_pixel[1] * weight;
                    data_pixel[2] += color_pixel[2] * weight;
                    data_pixel[3] += color_pixel[3] * weight;
                }
            }

            // Alpha hole mask is combination of other channels
            mask_pixel[3] += alpha_mask;

        }
    }

    // Fill Alpha Holes
    //
    // Alpha holes occur when an object becomes uncovered due to an object moving in front of it. Since
    // information is lost (the color of the pixels behind the front object) we have to make a guess at what
    // was there. But first we have to estimate where these holes are.
    //
    // As the pixels are blurred in the direction of movement, the original color alpha channel is sampled.
    // The estimate of the alpha *should be* is simply the fraction of the blurred pixel vector that is within
    // the mask of the non blurred color channel. As the blurred samples are summed for each pixel, the estimate will
    // converge to the desired alpha. You can then compare the desired alpha to the actual alpha to estimate the
    // holes.
    //
    // It's not intuitive at all, but it works. Turn on DEBUG_RENDER_ALPHA_MASK to see the result.
    //
    // Once we know what the alpha value of the hole (from the color alpha channel) is and the desired alpha (mask)
    // we can "boost" the color until the alpha equals the original value. NOTE: This technique cannot fix alphas
    // outside of the object but the artifacts are acceptable in motion.

#if !DEBUG_RENDER_ALPHA_MASK
    if (m_settings->fill_alpha_holes) {
#endif
        for (int y = 0; y < getHeight(); ++y) {
            for (int x = 0; x < getWidth(); ++x) {
                int index = INDEX_COL(x,y);
                float *color_pixel = color->getBuffer() + index;
                float *data_pixel = data + index;
                float *mask_pixel = mask + index;

                // The differences in alpha channels masked by the original alpha channel
                float fill_alpha = max(mask_pixel[3] - min(data_pixel[3], 1.0f),0.0f) * color_pixel[3];

#if DEBUG_RENDER_ALPHA_MASK
                data_pixel[0] = fill_alpha;
                data_pixel[1] = fill_alpha;
                data_pixel[2] = fill_alpha;
                data_pixel[3] = 1.0f;
#else
                if (fill_alpha <= 0.0F)
                    continue;

                // So if a color (color_pixel.rgb, alpha_mask) is the final color after the comp with holes, what
                // was the original color without the holes?
                //
                // if:      existing_alpha * f = existing_alpha + fill_alpha
                // then:    f = (existing_alpha + fill_alpha) / existing_alpha
                //
                // NOTE:    existing_alpha + fill_alpha != 1.0 which is the whole point of making the mask ;)

                float f = (data_pixel[3] + fill_alpha) / data_pixel[3];

                data_pixel[0] *= f;
                data_pixel[1] *= f;
                data_pixel[2] *= f;
                data_pixel[3] *= f;
#endif

            }
        }

#if !DEBUG_RENDER_ALPHA_MASK
    }
#endif

    // Delete mask
    MEM_freeN(mask);
}


MotionBlur2DOperation::DeepSample* MotionBlur2DOperation::alloc_sample (void)
{
    if (!samples) {
        // Allocate a bunch of samples
        DeepSampleBuffer *new_buffer = (DeepSampleBuffer *)MEM_callocN(sizeof(DeepSampleBuffer), "motion2D deep sample buffer");
        for (int i = 0; i < ARRAY_SIZE(new_buffer->buffers); ++i) {
            new_buffer->buffers[i].next_sample = samples;
            samples = &(new_buffer->buffers[i]);
        }

        // Save the buffer so we can delete it later
        new_buffer->next_buffer = buffers;
        buffers = new_buffer;
    }

    // Return the first sample
    DeepSample *new_sample = samples;
    samples = samples->next_sample;
    return new_sample;
}

bool MotionBlur2DOperation::deepSamplesSortFn(DeepSample *a, DeepSample *b) {
    return a->depth < b->depth;
}

void MotionBlur2DOperation::generateMotionBlurDeep(float *data, MemoryBuffer *color, MemoryBuffer *speed, MemoryBuffer *depth, MemoryBuffer *objid)
{
    // Fast sample allocators
    samples = NULL;
    buffers = NULL;

    // Allocate a per pixel sample buffer
    DeepSamplePixel *deep_samples = (DeepSamplePixel*) MEM_callocN(getWidth() * getHeight() * sizeof(DeepSamplePixel), "motion2D deep sample pixel buffer");

    // Clamp Multisample settings
    int multisample = m_settings->multisample;
    if (multisample < 0)        multisample = 1;
    else if (multisample > 8)   multisample = 8;
    int multisample2 = multisample * multisample;

    // Adjusted width and height given multisampling
    int width_multisample = getWidth() * multisample;
    int height_multisample = getHeight() * multisample;

    Sample line_samples[MAX_SAMPLES];
    int num_line_samples;

    // Samples Generation
    //
    // All this does is blurs the pixel contribution along the direction of movement for the pixel. If
    // multisampling is used then it acts internally like the images are X pixels larger (so
    // this gets slower) while taking care to keep the reading and writing of buffers correctly.
    //
    // Instead of mixing the samples directly as above, it records the samples into a linked list.
    //

    float forward_factor = m_settings->amount * multisample;
    float backwards_factor = m_settings->amount * multisample;

    for (int y = 0; y < height_multisample; ++y) {
        for (int x = 0; x < width_multisample; ++x) {
            int xm = x/multisample;
            int ym = y/multisample;

            // Record sample
            int index_col = INDEX_COL(xm,ym);
            float *color_pixel = color->getBuffer() + index_col;

            // Skip full transparent pixels
            if (color_pixel[3] == 0.0f)
                continue;

            int index_vec = INDEX_COL(xm,ym);
            int index_val = INDEX_VAL(xm,ym);
            float *speed_pixel = speed->getBuffer() + index_vec;
            float *depth_pixel = depth->getBuffer() + index_val;
            float *objid_pixel = objid->getBuffer() + index_val;

            // Calculate a motion vector
            float z = *depth_pixel;

            int x0 = x + speed_pixel[0] * forward_factor;
            int y0 = y + speed_pixel[1] * forward_factor;
            int x1 = x - speed_pixel[0] * backwards_factor;
            int y1 = y - speed_pixel[1] * backwards_factor;

            // Build blur samples
            line (x0, y0, x1, y1, line_samples, &num_line_samples);
            float weight = 1.0f / (num_line_samples * multisample2);

            for (int s = 0; s < num_line_samples; ++s) {
                int xs = line_samples[s].x;
                int ys = line_samples[s].y;

                // Alpha ramp
                float alpha = (float) (s+1) / (float) num_line_samples;
                if (num_line_samples > 1)
                    alpha = (alpha < 0.5F) ? (alpha*2.0F) : (1.0F - (alpha-0.5F)*2.0F);    //a=0,r=0; a=0.5,r=1.0; a=1.0,r=0.0
                alpha *= color_pixel[3];

                // If outside image, ignore
                if (xs < 0 || xs >= width_multisample || ys < 0 || ys >= height_multisample) {
                    // Do nothing

                } else {
                    int xm_s = xs/multisample;
                    int ym_s = ys/multisample;

                    // Record sample
                    DeepSamplePixel &pixel_samples = *(deep_samples + INDEX_VAL(xm_s,ym_s));

                    DeepSample *sample = pixel_samples.samples;
                    DeepSample *prev_sample = NULL;

                    // Search for sample with same ID
                    while (sample) {
                        if (sample->obj_id == *objid_pixel) {

                            // Move sample to front of list if not already (like a MRU cache)
                            if (prev_sample) {
                                // Remove sample
                                prev_sample->next_sample = sample->next_sample;

                                // Push to front
                                sample->next_sample = pixel_samples.samples;
                                pixel_samples.samples = sample;
                            }

                            break;
                        }

                        prev_sample=sample;
                        sample=sample->next_sample;
                    }

                    // Build or update the sample
                    if (!sample) {
                        sample = alloc_sample();
                        sample->obj_id = *objid_pixel;
                        sample->color[0] = color_pixel[0] * weight;
                        sample->color[1] = color_pixel[1] * weight;
                        sample->color[2] = color_pixel[2] * weight;
                        sample->color[3] = color_pixel[3] * weight;
                        sample->depth = z;
                        sample->max_alpha = alpha;

                        sample->next_sample = pixel_samples.samples;
                        pixel_samples.samples = sample;
                    } else {
                        sample->color[0] += color_pixel[0] * weight;
                        sample->color[1] += color_pixel[1] * weight;
                        sample->color[2] += color_pixel[2] * weight;
                        sample->color[3] += color_pixel[3] * weight;
                        sample->depth = std::min(sample->depth, z);
                        sample->max_alpha = std::max(sample->max_alpha, alpha);
                    }

                }
            }
        }
    }

    //
    // Compositing of samples
    //

    // Copy samples into array for easier access
    std::vector<DeepSample*> samples_array;

    for (int y = 0; y < getHeight(); ++y) {
        for (int x = 0; x < getWidth(); ++x) {

            DeepSamplePixel &pixel_samples = *(deep_samples + INDEX_VAL(x, y));
            float *data_pixel = data + INDEX_COL(x,y);

            // Copy to array
            samples_array.clear();
            DeepSample *sample = pixel_samples.samples;
            while (sample) {
                samples_array.push_back(sample);
                sample = sample->next_sample;
            }

            // Sort samples based on depth. Front samples first.
            std::sort(samples_array.begin(), samples_array.end(), MotionBlur2DOperation::deepSamplesSortFn);

            // Combine samples from same object
            for (int i = 0; i < samples_array.size(); ++i) {
                sample = samples_array[i];

                // Correct alpha to calculated alpha
                if (sample->color[3] <= 0.0F || sample->max_alpha <= 0.0F)
                    continue;

                // Normalize sample and premultiply alphas
                sample->color[0] = sample->color[0] / sample->color[3] * sample->max_alpha;
                sample->color[1] = sample->color[1] / sample->color[3] * sample->max_alpha;
                sample->color[2] = sample->color[2] / sample->color[3] * sample->max_alpha;
                sample->color[3] = sample->max_alpha;

                // Comp samples behind existing samples so far. This is an "over" operation.
                float front_mix = 1.0f;
                float back_mix = (1.0f - data_pixel[3]);

                data_pixel[0] = data_pixel[0] * front_mix + sample->color[0] * back_mix;
                data_pixel[1] = data_pixel[1] * front_mix + sample->color[1] * back_mix;
                data_pixel[2] = data_pixel[2] * front_mix + sample->color[2] * back_mix;
                data_pixel[3] = data_pixel[3] * front_mix + sample->color[3] * back_mix;
            }

            // Un-premultiply Alphas
//            data_pixel[0] = data_pixel[0] / data_pixel[3];
//            data_pixel[1] = data_pixel[1] / data_pixel[3];
//            data_pixel[2] = data_pixel[2] / data_pixel[3];
        }
    }

    // Delete all sample buffers
    while (buffers) {
        DeepSampleBuffer *next_buffer = buffers->next_buffer;
        MEM_freeN(buffers);
        buffers = next_buffer;
    }

    MEM_freeN(deep_samples);
}



