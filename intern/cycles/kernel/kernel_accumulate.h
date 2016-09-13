/*
 * Copyright 2011-2013 Blender Foundation
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

CCL_NAMESPACE_BEGIN

/* BSDF Eval
 *
 * BSDF evaluation result, split per BSDF type. This is used to accumulate
 * render passes separately. */

ccl_device_inline void bsdf_eval_init(BsdfEval *eval, ClosureType type, float3 value, int use_light_pass)
{
#ifdef __PASSES__
	eval->use_light_pass = use_light_pass;

	if(eval->use_light_pass) {
		eval->diffuse = make_float3(0.0f, 0.0f, 0.0f);
		eval->glossy = make_float3(0.0f, 0.0f, 0.0f);
		eval->transmission = make_float3(0.0f, 0.0f, 0.0f);
		eval->transparent = make_float3(0.0f, 0.0f, 0.0f);
		eval->subsurface = make_float3(0.0f, 0.0f, 0.0f);
		eval->scatter = make_float3(0.0f, 0.0f, 0.0f);

		if(type == CLOSURE_BSDF_TRANSPARENT_ID)
			eval->transparent = value;
		else if(CLOSURE_IS_BSDF_DIFFUSE(type))
			eval->diffuse = value;
		else if(CLOSURE_IS_BSDF_GLOSSY(type))
			eval->glossy = value;
		else if(CLOSURE_IS_BSDF_TRANSMISSION(type))
			eval->transmission = value;
		else if(CLOSURE_IS_BSDF_BSSRDF(type))
			eval->subsurface = value;
		else if(CLOSURE_IS_PHASE(type))
			eval->scatter = value;
	}
	else
		eval->diffuse = value;
#else
	*eval = value;
#endif
}

/* TODO(sergey): This is just a workaround for annoying 6.5 compiler bug. */
#if !defined(__KERNEL_CUDA__) || __CUDA_ARCH__ < 500
ccl_device_inline
#else
ccl_device_noinline
#endif
void bsdf_eval_accum(BsdfEval *eval, ClosureType type, float3 value)
{
#ifdef __PASSES__
	if(eval->use_light_pass) {
		if(CLOSURE_IS_BSDF_DIFFUSE(type))
			eval->diffuse += value;
		else if(CLOSURE_IS_BSDF_GLOSSY(type))
			eval->glossy += value;
		else if(CLOSURE_IS_BSDF_TRANSMISSION(type))
			eval->transmission += value;
		else if(CLOSURE_IS_BSDF_BSSRDF(type))
			eval->subsurface += value;
		else if(CLOSURE_IS_PHASE(type))
			eval->scatter += value;

		/* skipping transparent, this function is used by for eval(), will be zero then */
	}
	else
		eval->diffuse += value;
#else
	*eval += value;
#endif
}

ccl_device_inline bool bsdf_eval_is_zero(BsdfEval *eval)
{
#ifdef __PASSES__
	if(eval->use_light_pass) {
		return is_zero(eval->diffuse)
			&& is_zero(eval->glossy)
			&& is_zero(eval->transmission)
			&& is_zero(eval->transparent)
			&& is_zero(eval->subsurface)
			&& is_zero(eval->scatter);
	}
	else
		return is_zero(eval->diffuse);
#else
	return is_zero(*eval);
#endif
}

ccl_device_inline void bsdf_eval_mul(BsdfEval *eval, float3 value)
{
#ifdef __PASSES__
	if(eval->use_light_pass) {
		eval->diffuse *= value;
		eval->glossy *= value;
		eval->transmission *= value;
		eval->subsurface *= value;
		eval->scatter *= value;

		/* skipping transparent, this function is used by for eval(), will be zero then */
	}
	else
		eval->diffuse *= value;
#else
	*eval *= value;
#endif
}

/* Path Radiance
 *
 * We accumulate different render passes separately. After summing at the end
 * to get the combined result, it should be identical. We definite directly
 * visible as the first non-transparent hit, while indirectly visible are the
 * bounces after that. */

ccl_device_inline void path_radiance_init(PathRadiance *L, int use_light_pass)
{
#ifdef __PASSES__
    L->use_light_pass = use_light_pass;
    L->ao = make_float3(0.0f, 0.0f, 0.0f);
    L->mist = 0.0f;

	/* clear all */
    for (int i = 0; i < MAX_NUM_LIGHT_BUFFERS; ++i) {
        if(use_light_pass) {
            PathRadianceBuffer *b = &(L->buffers[i]);

            b->background = make_float3(0.0f, 0.0f, 0.0f);

            b->emission = make_float3(0.0f, 0.0f, 0.0f);
            b->indirect = make_float3(0.0f, 0.0f, 0.0f);
            b->direct_throughput = make_float3(0.0f, 0.0f, 0.0f);
            b->direct_emission = make_float3(0.0f, 0.0f, 0.0f);

            b->color_diffuse = make_float3(0.0f, 0.0f, 0.0f);
            b->color_glossy = make_float3(0.0f, 0.0f, 0.0f);
            b->color_transmission = make_float3(0.0f, 0.0f, 0.0f);
            b->color_subsurface = make_float3(0.0f, 0.0f, 0.0f);
            b->color_scatter = make_float3(0.0f, 0.0f, 0.0f);

            b->direct_diffuse = make_float3(0.0f, 0.0f, 0.0f);
            b->direct_glossy = make_float3(0.0f, 0.0f, 0.0f);
            b->direct_transmission = make_float3(0.0f, 0.0f, 0.0f);
            b->direct_subsurface = make_float3(0.0f, 0.0f, 0.0f);
            b->direct_scatter = make_float3(0.0f, 0.0f, 0.0f);

            b->indirect_diffuse = make_float3(0.0f, 0.0f, 0.0f);
            b->indirect_glossy = make_float3(0.0f, 0.0f, 0.0f);
            b->indirect_transmission = make_float3(0.0f, 0.0f, 0.0f);
            b->indirect_subsurface = make_float3(0.0f, 0.0f, 0.0f);
            b->indirect_scatter = make_float3(0.0f, 0.0f, 0.0f);

            b->path_diffuse = make_float3(0.0f, 0.0f, 0.0f);
            b->path_glossy = make_float3(0.0f, 0.0f, 0.0f);
            b->path_transmission = make_float3(0.0f, 0.0f, 0.0f);
            b->path_subsurface = make_float3(0.0f, 0.0f, 0.0f);
            b->path_scatter = make_float3(0.0f, 0.0f, 0.0f);

            b->shadow = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
        } else {
            PathRadianceBuffer *b = &(L->buffers[i]);
            b->emission = make_float3(0.0f, 0.0f, 0.0f);
        }
    }
#else
    *L = make_float3(0.0f, 0.0f, 0.0f);
#endif

}

ccl_device_inline void path_radiance_init(PathRadianceAccum *L, int use_light_pass)
{
	/* clear all */
#ifdef __PASSES__
    L->use_light_pass = use_light_pass;
    L->background = make_float3(0.0f, 0.0f, 0.0f);
    L->ao = make_float3(0.0f, 0.0f, 0.0f);
    L->mist = 0.0f;

    if(use_light_pass) {
        L->emission = make_float3(0.0f, 0.0f, 0.0f);
        L->indirect = make_float3(0.0f, 0.0f, 0.0f);
        L->direct_throughput = make_float3(0.0f, 0.0f, 0.0f);
        L->direct_emission = make_float3(0.0f, 0.0f, 0.0f);

        L->color_diffuse = make_float3(0.0f, 0.0f, 0.0f);
        L->color_glossy = make_float3(0.0f, 0.0f, 0.0f);
        L->color_transmission = make_float3(0.0f, 0.0f, 0.0f);
        L->color_subsurface = make_float3(0.0f, 0.0f, 0.0f);
        L->color_scatter = make_float3(0.0f, 0.0f, 0.0f);

        L->direct_diffuse = make_float3(0.0f, 0.0f, 0.0f);
        L->direct_glossy = make_float3(0.0f, 0.0f, 0.0f);
        L->direct_transmission = make_float3(0.0f, 0.0f, 0.0f);
        L->direct_subsurface = make_float3(0.0f, 0.0f, 0.0f);
        L->direct_scatter = make_float3(0.0f, 0.0f, 0.0f);

        L->indirect_diffuse = make_float3(0.0f, 0.0f, 0.0f);
        L->indirect_glossy = make_float3(0.0f, 0.0f, 0.0f);
        L->indirect_transmission = make_float3(0.0f, 0.0f, 0.0f);
        L->indirect_subsurface = make_float3(0.0f, 0.0f, 0.0f);
        L->indirect_scatter = make_float3(0.0f, 0.0f, 0.0f);

        L->path_diffuse = make_float3(0.0f, 0.0f, 0.0f);
        L->path_glossy = make_float3(0.0f, 0.0f, 0.0f);
        L->path_transmission = make_float3(0.0f, 0.0f, 0.0f);
        L->path_subsurface = make_float3(0.0f, 0.0f, 0.0f);
        L->path_scatter = make_float3(0.0f, 0.0f, 0.0f);

        L->shadow = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
    } else {
        L->emission = make_float3(0.0f, 0.0f, 0.0f);
    }
#else
    *L = make_float3(0.0f, 0.0f, 0.0f);
#endif
}

ccl_device_inline void path_radiance_accum_buffers(PathRadiance *L, PathRadianceAccum *Laccum)
{
    Laccum->ao += L->ao;
    Laccum->mist += L->mist;

    for (int i = 0; i < MAX_NUM_LIGHT_BUFFERS; ++i) {
        Laccum->background += L->buffers[i].background;

        Laccum->emission += L->buffers[i].emission;
        Laccum->indirect += L->buffers[i].indirect;
        Laccum->direct_throughput += L->buffers[i].direct_throughput;
        Laccum->direct_emission += L->buffers[i].direct_emission;

        Laccum->color_diffuse += L->buffers[i].color_diffuse;
        Laccum->color_glossy += L->buffers[i].color_glossy;
        Laccum->color_transmission += L->buffers[i].color_transmission;
        Laccum->color_subsurface += L->buffers[i].color_subsurface;
        Laccum->color_scatter += L->buffers[i].color_scatter;

        Laccum->direct_diffuse += L->buffers[i].direct_diffuse;
        Laccum->direct_glossy += L->buffers[i].direct_glossy;
        Laccum->direct_transmission += L->buffers[i].direct_transmission;
        Laccum->direct_subsurface += L->buffers[i].direct_subsurface;
        Laccum->direct_scatter += L->buffers[i].direct_scatter;

        Laccum->indirect_diffuse += L->buffers[i].indirect_diffuse;
        Laccum->indirect_glossy += L->buffers[i].indirect_glossy;
        Laccum->indirect_transmission += L->buffers[i].indirect_transmission;
        Laccum->indirect_subsurface += L->buffers[i].indirect_subsurface;
        Laccum->indirect_scatter += L->buffers[i].indirect_scatter;

        Laccum->path_diffuse += L->buffers[i].path_diffuse;
        Laccum->path_glossy += L->buffers[i].path_glossy;
        Laccum->path_transmission += L->buffers[i].path_transmission;
        Laccum->path_subsurface += L->buffers[i].path_subsurface;
        Laccum->path_scatter += L->buffers[i].path_scatter;

        Laccum->shadow += L->buffers[i].shadow;
    }

}

ccl_device_inline void path_radiance_bsdf_bounce(PathRadiance *L, int nL, ccl_addr_space float3 *throughput,
	BsdfEval *bsdf_eval, float bsdf_pdf, int bounce, int bsdf_label)
{
	float inverse_pdf = 1.0f/bsdf_pdf;

#ifdef __PASSES__
	if(L->use_light_pass) {
		if(bounce == 0 && !(bsdf_label & LABEL_TRANSPARENT)) {
			/* first on directly visible surface */
			float3 value = *throughput*inverse_pdf;
            PathRadianceBuffer *b = &(L->buffers[nL]);

			b->path_diffuse = bsdf_eval->diffuse*value;
			b->path_glossy = bsdf_eval->glossy*value;
			b->path_transmission = bsdf_eval->transmission*value;
			b->path_subsurface = bsdf_eval->subsurface*value;
			b->path_scatter = bsdf_eval->scatter*value;

			*throughput = b->path_diffuse + b->path_glossy + b->path_transmission + b->path_subsurface + b->path_scatter;
			
			b->direct_throughput = *throughput;
		}
		else {
			/* transparent bounce before first hit, or indirectly visible through BSDF */
			float3 sum = (bsdf_eval->diffuse + bsdf_eval->glossy + bsdf_eval->transmission + bsdf_eval->transparent +
						  bsdf_eval->subsurface + bsdf_eval->scatter) * inverse_pdf;
			*throughput *= sum;
		}
	}
	else
		*throughput *= bsdf_eval->diffuse*inverse_pdf;
#else
	*throughput *= *bsdf_eval*inverse_pdf;
#endif
}

ccl_device_inline void path_radiance_accum_emission(PathRadiance *L, float3 throughput, float3 emission[MAX_NUM_LIGHT_BUFFERS], int bounce)
{
    for (int i = 0; i < MAX_NUM_LIGHT_BUFFERS; ++i) {
#ifdef __PASSES__
        if(L->use_light_pass) {
            if(bounce == 0)
                L->buffers[i].emission += throughput*emission[i];
            else if(bounce == 1)
                L->buffers[i].direct_emission += throughput*emission[i];
            else
                L->buffers[i].indirect += throughput*emission[i];
        }
        else
            L->buffers[i].emission += throughput*emission[i];
#else
        L->buffers[i] += throughput*emission[i];
#endif
    }
}

ccl_device_inline void path_radiance_accum_emission(PathRadiance *L, int nL, float3 throughput, float3 emission, int bounce)
{
#ifdef __PASSES__
    if(L->use_light_pass) {
        if(bounce == 0)
            L->buffers[nL].emission += throughput*emission;
        else if(bounce == 1)
            L->buffers[nL].direct_emission += throughput*emission;
        else
            L->buffers[nL].indirect += throughput*emission;
    }
    else
        L->buffers[nL].emission += throughput*emission;
#else
    *L += throughput*emission;
#endif
}

ccl_device_inline void path_radiance_accum_ao(PathRadiance *L, int nL, float3 throughput, float3 alpha, float3 bsdf, float3 ao, int bounce)
{
#ifdef __PASSES__
	if(L->use_light_pass) {
		if(bounce == 0) {
			/* directly visible lighting */
			L->buffers[nL].direct_diffuse += throughput*bsdf*ao;
			L->ao += alpha*throughput*ao;
		}
		else {
			/* indirectly visible lighting after BSDF bounce */
			L->buffers[nL].indirect += throughput*bsdf*ao;
		}
	}
	else
		L->buffers[nL].emission += throughput*bsdf*ao;
#else
	L->buffers[nL] += throughput*bsdf*ao;
#endif
}

ccl_device_inline void path_radiance_accum_light(PathRadiance *L, int nL, float3 throughput, BsdfEval *bsdf_eval, float3 shadow, float shadow_fac, int bounce, bool is_lamp)
{
#ifdef __PASSES__
	if(L->use_light_pass) {
		if(bounce == 0) {
			/* directly visible lighting */
			L->buffers[nL].direct_diffuse += throughput*bsdf_eval->diffuse*shadow;
			L->buffers[nL].direct_glossy += throughput*bsdf_eval->glossy*shadow;
			L->buffers[nL].direct_transmission += throughput*bsdf_eval->transmission*shadow;
			L->buffers[nL].direct_subsurface += throughput*bsdf_eval->subsurface*shadow;
			L->buffers[nL].direct_scatter += throughput*bsdf_eval->scatter*shadow;

			if(is_lamp) {
				L->buffers[nL].shadow.x += shadow.x*shadow_fac;
				L->buffers[nL].shadow.y += shadow.y*shadow_fac;
				L->buffers[nL].shadow.z += shadow.z*shadow_fac;
			}
		}
		else {
			/* indirectly visible lighting after BSDF bounce */
			float3 sum = bsdf_eval->diffuse + bsdf_eval->glossy + bsdf_eval->transmission + bsdf_eval->subsurface + bsdf_eval->scatter;
			L->buffers[nL].indirect += throughput*sum*shadow;
		}
	}
	else
		L->buffers[nL].emission += throughput*bsdf_eval->diffuse*shadow;
#else
	L[nL] += throughput*(*bsdf_eval)*shadow;
#endif
}

ccl_device_inline void path_radiance_accum_background(PathRadiance *L, int nL, float3 throughput, float3 value, int bounce)
{
#ifdef __PASSES__
	if(L->use_light_pass) {
		if(bounce == 0)
			L->buffers[nL].background += throughput*value;
		else if(bounce == 1)
			L->buffers[nL].direct_emission += throughput*value;
		else
			L->buffers[nL].indirect += throughput*value;
	}
	else
		L->buffers[nL].emission += throughput*value;
#else
	L[nL] += throughput*value;
#endif
}

ccl_device_inline void path_radiance_sum_indirect(PathRadiance *L, int nL)
{
#ifdef __PASSES__
	/* this division is a bit ugly, but means we only have to keep track of
	 * only a single throughput further along the path, here we recover just
	 * the indirect path that is not influenced by any particular BSDF type */
	if(L->use_light_pass) {
        PathRadianceBuffer *b = &(L->buffers[nL]);

		b->direct_emission = safe_divide_color(b->direct_emission, b->direct_throughput);
		b->direct_diffuse += b->path_diffuse*b->direct_emission;
		b->direct_glossy += b->path_glossy*b->direct_emission;
		b->direct_transmission += b->path_transmission*b->direct_emission;
		b->direct_subsurface += b->path_subsurface*b->direct_emission;
		b->direct_scatter += b->path_scatter*b->direct_emission;

		b->indirect = safe_divide_color(b->indirect, b->direct_throughput);
		b->indirect_diffuse += b->path_diffuse*b->indirect;
		b->indirect_glossy += b->path_glossy*b->indirect;
		b->indirect_transmission += b->path_transmission*b->indirect;
		b->indirect_subsurface += b->path_subsurface*b->indirect;
		b->indirect_scatter += b->path_scatter*b->indirect;
	}
#endif
}

ccl_device_inline void path_radiance_sum_indirect(PathRadianceAccum *L)
{
#ifdef __PASSES__
	/* this division is a bit ugly, but means we only have to keep track of
	 * only a single throughput further along the path, here we recover just
	 * the indirect path that is not influenced by any particular BSDF type */
	if(L->use_light_pass) {
		L->direct_emission = safe_divide_color(L->direct_emission, L->direct_throughput);
		L->direct_diffuse += L->path_diffuse*L->direct_emission;
		L->direct_glossy += L->path_glossy*L->direct_emission;
		L->direct_transmission += L->path_transmission*L->direct_emission;
		L->direct_subsurface += L->path_subsurface*L->direct_emission;
		L->direct_scatter += L->path_scatter*L->direct_emission;

		L->indirect = safe_divide_color(L->indirect, L->direct_throughput);
		L->indirect_diffuse += L->path_diffuse*L->indirect;
		L->indirect_glossy += L->path_glossy*L->indirect;
		L->indirect_transmission += L->path_transmission*L->indirect;
		L->indirect_subsurface += L->path_subsurface*L->indirect;
		L->indirect_scatter += L->path_scatter*L->indirect;
	}
#endif
}

ccl_device_inline void path_radiance_reset_indirect(PathRadiance *L, int nL)
{
#ifdef __PASSES__
	if(L->use_light_pass) {
        PathRadianceBuffer *b = &(L->buffers[nL]);

		b->path_diffuse = make_float3(0.0f, 0.0f, 0.0f);
		b->path_glossy = make_float3(0.0f, 0.0f, 0.0f);
		b->path_transmission = make_float3(0.0f, 0.0f, 0.0f);
		b->path_subsurface = make_float3(0.0f, 0.0f, 0.0f);
		b->path_scatter = make_float3(0.0f, 0.0f, 0.0f);

		b->direct_emission = make_float3(0.0f, 0.0f, 0.0f);
		b->indirect = make_float3(0.0f, 0.0f, 0.0f);
	}
#endif
}

ccl_device_inline float3 path_radiance_clamp_and_sum(KernelGlobals *kg, PathRadianceAccum *L)
{
	float3 L_sum;
	/* Light Passes are used */
#ifdef __PASSES__
	float3 L_direct, L_indirect;
	float clamp_direct = kernel_data.integrator.sample_clamp_direct;
	float clamp_indirect = kernel_data.integrator.sample_clamp_indirect;
	if(L->use_light_pass) {
		path_radiance_sum_indirect(L);

		L_direct = L->direct_diffuse + L->direct_glossy + L->direct_transmission + L->direct_subsurface + L->direct_scatter + L->emission;
		L_indirect = L->indirect_diffuse + L->indirect_glossy + L->indirect_transmission + L->indirect_subsurface + L->indirect_scatter;

		if(!kernel_data.background.transparent)
			L_direct += L->background;

		L_sum = L_direct + L_indirect;
		float sum = fabsf((L_sum).x) + fabsf((L_sum).y) + fabsf((L_sum).z);

		/* Reject invalid value */
		if(!isfinite(sum)) {
			L_sum = make_float3(0.0f, 0.0f, 0.0f);

			L->direct_diffuse = make_float3(0.0f, 0.0f, 0.0f);
			L->direct_glossy = make_float3(0.0f, 0.0f, 0.0f);
			L->direct_transmission = make_float3(0.0f, 0.0f, 0.0f);
			L->direct_subsurface = make_float3(0.0f, 0.0f, 0.0f);
			L->direct_scatter = make_float3(0.0f, 0.0f, 0.0f);

			L->indirect_diffuse = make_float3(0.0f, 0.0f, 0.0f);
			L->indirect_glossy = make_float3(0.0f, 0.0f, 0.0f);
			L->indirect_transmission = make_float3(0.0f, 0.0f, 0.0f);
			L->indirect_subsurface = make_float3(0.0f, 0.0f, 0.0f);
			L->indirect_scatter = make_float3(0.0f, 0.0f, 0.0f);

			L->emission = make_float3(0.0f, 0.0f, 0.0f);
		}

		/* Clamp direct and indirect samples */
#ifdef __CLAMP_SAMPLE__
		else if(sum > clamp_direct || sum > clamp_indirect) {
			float scale;

			/* Direct */
			float sum_direct = fabsf(L_direct.x) + fabsf(L_direct.y) + fabsf(L_direct.z);
			if(sum_direct > clamp_direct) {
				scale = clamp_direct/sum_direct;
				L_direct *= scale;

				L->direct_diffuse *= scale;
				L->direct_glossy *= scale;
				L->direct_transmission *= scale;
				L->direct_subsurface *= scale;
				L->direct_scatter *= scale;
				L->emission *= scale;
				L->background *= scale;
			}

			/* Indirect */
			float sum_indirect = fabsf(L_indirect.x) + fabsf(L_indirect.y) + fabsf(L_indirect.z);
			if(sum_indirect > clamp_indirect) {
				scale = clamp_indirect/sum_indirect;
				L_indirect *= scale;

				L->indirect_diffuse *= scale;
				L->indirect_glossy *= scale;
				L->indirect_transmission *= scale;
				L->indirect_subsurface *= scale;
				L->indirect_scatter *= scale;
			}

			/* Sum again, after clamping */
			L_sum = L_direct + L_indirect;
		}
#endif

		return L_sum;
	}

	/* No Light Passes */
	else
		L_sum = L->emission;
#else
	L_sum = *L;
#endif

	/* Reject invalid value */
	float sum = fabsf((L_sum).x) + fabsf((L_sum).y) + fabsf((L_sum).z);
	if(!isfinite(sum))
		L_sum = make_float3(0.0f, 0.0f, 0.0f);

	return L_sum;
}

ccl_device_inline void path_radiance_accum_sample(PathRadiance *L, PathRadiance *L_sample, int nL, int num_samples)
{
	float fac = 1.0f/num_samples;

#ifdef __PASSES__
    PathRadianceBuffer *b = &(L->buffers[nL]);
    PathRadianceBuffer *b_sample = &(L_sample->buffers[nL]);

	b->direct_diffuse += b_sample->direct_diffuse*fac;
	b->direct_glossy += b_sample->direct_glossy*fac;
	b->direct_transmission += b_sample->direct_transmission*fac;
	b->direct_subsurface += b_sample->direct_subsurface*fac;
	b->direct_scatter += b_sample->direct_scatter*fac;

	b->indirect_diffuse += b_sample->indirect_diffuse*fac;
	b->indirect_glossy += b_sample->indirect_glossy*fac;
	b->indirect_transmission += b_sample->indirect_transmission*fac;
	b->indirect_subsurface += b_sample->indirect_subsurface*fac;
	b->indirect_scatter += b_sample->indirect_scatter*fac;

	b->emission += b_sample->emission*fac;
	b->shadow += b_sample->shadow*fac;
	b->background += b_sample->background*fac;

	L->ao += L_sample->ao*fac;
	L->mist += L_sample->mist*fac;
#else
	*L += *L_sample * fac;
#endif
}

CCL_NAMESPACE_END

