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

#define INDEX(x,y) ((y * this->getWidth() + x) * COM_NUM_CHANNELS_COLOR)
#define MAX_SAMPLES 1024*16
#define DEBUG_RENDER_ALPHA_MASK 0

MotionBlur2DOperation::MotionBlur2DOperation() : NodeOperation()
{
	this->addInputSocket(COM_DT_COLOR);
	this->addInputSocket(COM_DT_COLOR);
	this->addOutputSocket(COM_DT_COLOR);
	this->m_settings = NULL;
	this->m_inputImageProgram = NULL;
	this->m_inputSpeedProgram = NULL;
	this->m_cachedInstance = NULL;
	setComplex(true);
}
void MotionBlur2DOperation::initExecution()
{
	initMutex();
	this->m_inputImageProgram = getInputSocketReader(0);
	this->m_inputSpeedProgram = getInputSocketReader(1);
	this->m_cachedInstance = NULL;
}

void MotionBlur2DOperation::deinitExecution()
{
	deinitMutex();
	if (this->m_cachedInstance) {
		MEM_freeN(this->m_cachedInstance);
		this->m_cachedInstance = NULL;
	}
}

void *MotionBlur2DOperation::initializeTileData(rcti *rect)
{
	if (this->m_cachedInstance) {
		return this->m_cachedInstance;
	}
	
	lockMutex();
	if (this->m_cachedInstance == NULL) {
		MemoryBuffer *color = (MemoryBuffer *)this->m_inputImageProgram->initializeTileData(rect);
		MemoryBuffer *speed = (MemoryBuffer *)this->m_inputSpeedProgram->initializeTileData(rect);
		float *data = (float *)MEM_callocN(MEM_allocN_len(color->getBuffer()), "motion2D data buffer");
		this->generateMotionBlur(data, color, speed);
		this->m_cachedInstance = data;
	}
	unlockMutex();
	return this->m_cachedInstance;
}

void MotionBlur2DOperation::executePixel(float output[4], int x, int y, void *data)
{
	float *buffer = (float *) data;
	copy_v4_v4(output, &buffer[INDEX(x,y)]);
}

bool MotionBlur2DOperation::determineDependingAreaOfInterest(rcti *input, ReadBufferOperation *readOperation, rcti *output)
{
	if (this->m_cachedInstance == NULL) {
		rcti newInput;
		newInput.xmax = this->getWidth();
		newInput.xmin = 0;
		newInput.ymax = this->getHeight();
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
    int multisample = this->m_settings->multisample;
    if (multisample < 0)        multisample = 1;
    else if (multisample > 4)   multisample = 4;
    int multisample2 = multisample * multisample;

    // Adjusted width and height given multisampling
    int width_multisample = this->getWidth() * multisample;
    int height_multisample = this->getHeight() * multisample;

    Sample samples[MAX_SAMPLES];
    int num_samples;

    // Blurring of pixels
    //
    // All this does is blurs the pixel contribution along the direction of movement for the pixel. If
    // multisampling is used then it acts internally like the images are X pixels larger (so
    // this gets slower) while taking care to keep the reading and writing of buffers correctly.
    //

    float forward_factor = (this->m_settings->blur_forwards * this->m_settings->amount) * multisample;
    float backwards_factor = (this->m_settings->blur_backwards * this->m_settings->amount) * multisample;

    for (int y = 0; y < height_multisample; ++y) {
        for (int x = 0; x < width_multisample; ++x) {
            // Render sample
            int index = INDEX(x/multisample,y/multisample);
            float *color_pixel = color->getBuffer() + index;

            // Skip full transparent pixels
            if (color_pixel[3] == 0.0f)
                continue;

            float *speed_pixel = speed->getBuffer() + index;
            float *mask_pixel = mask + index;

            // Calculate a motion vector

            int x0 = x + speed_pixel[0] * forward_factor;
            int y0 = y + speed_pixel[1] * forward_factor;
            int x1 = x - speed_pixel[0] * backwards_factor;
            int y1 = y - speed_pixel[1] * backwards_factor;

            // Build blur samples
            line (x0, y0, x1, y1, samples, &num_samples);
            float weight = 1.0f / (num_samples * multisample2);

            // Alpha hole mask
            float alpha_mask = 0.0f;

            for (int s = 0; s < num_samples; ++s) {
                int xs = samples[s].x;
                int ys = samples[s].y;

                // If outside image, assume alpha = 1
                if (xs < 0 || xs >= width_multisample || ys < 0 || ys >= height_multisample) {
                    alpha_mask += weight;

                } else {
                    int index_s = INDEX(xs/multisample,ys/multisample);

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
    if (this->m_settings->fill_alpha_holes) {
#endif
        for (int y = 0; y < getHeight(); ++y) {
            for (int x = 0; x < getWidth(); ++x) {
                int index = INDEX(x,y);
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



