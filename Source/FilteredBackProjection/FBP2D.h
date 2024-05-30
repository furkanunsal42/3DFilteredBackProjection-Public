#pragma once

#include <memory>

#include "ComputeProgram.h"
#include "Texture2D.h"

class FBP2D {
public:
	enum FilterType {
		RAM_LAK		= 0,
		SHEPP_LOGAN = 1,
		COSINE		= 2,
	};

	// main functions
	void project_forward(Texture2D& source_image, Texture2D& target_image, float total_rotation_count);
	void project_forward_equiangular(Texture2D& source_image, Texture2D& target_image, float xray_source_distance, float detector_panel_distance, float detector_width, float volume_width, float total_rotation_count);
	void project_backward(Texture2D& sinogram, Texture2D& out_texture);
	
	std::shared_ptr<Texture2D> project_forward(Texture2D& source_image, int projection_count, float total_rotation_count);
	std::shared_ptr<Texture2D> project_forward_equiangular(Texture2D& source_image, int projection_count, float xray_source_distance, float detector_panel_distance, float detector_width, float volume_width, float total_rotation_count);
	std::shared_ptr<Texture2D> project_backward(Texture2D& sinogram);
	
	void parallelize_fan_beam_sinogram(Texture2D& source_sinogram, Texture2D& target_sinogram, float xray_source_distance, float detector_panel_distance, float detector_width, float volume_width, float total_rotation_count);

	void apply_filter(Texture2D& sinogram, FilterType filter_type = FilterType::RAM_LAK);

	// utility functions
	void _texture_blit_int32_to_float(Texture2D& source_int_texture, Texture2D& target_float_texture);
	void _texture_blit_float_to_int32(Texture2D& source_float_texture, Texture2D& target_int_texture);
	void _texture_blit_float1_to_float4(Texture2D& source_float_texture, Texture2D& target_float_texture);
	void _texture_blit_float1_to_float1(Texture2D& source_float_texture, Texture2D& target_float_texture);
	void _texture_blit_float1_to_r16(Texture2D& source_float_texture, Texture2D& target_float_texture);


	void _int_texture_clip_overflow(Texture2D& int_texture);

	void _int_texture_add();
	void _int_texture_mult();
	void _int_texture_div(Texture2D& int_texture, float denominator);

	// compute shaders
	std::shared_ptr<ComputeProgram> cp_radon_transform_int				= std::make_shared<ComputeProgram>(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/radon_transform_int.comp"));
	std::shared_ptr<ComputeProgram> cp_radon_transform_float			= std::make_shared<ComputeProgram>(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/radon_transform_float.comp"));
	std::shared_ptr<ComputeProgram> cp_radon_transform_inverse_int		= std::make_shared<ComputeProgram>(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/radon_transform_inverse_int.comp"));
	std::shared_ptr<ComputeProgram> cp_radon_transform_inverse_float	= std::make_shared<ComputeProgram>(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/radon_transform_inverse_float.comp"));

	std::shared_ptr<ComputeProgram> cp_radon_transform_fan_equiangular_float			 = std::make_shared<ComputeProgram>(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/radon_transform_fan_equiangular_float.comp"));
	std::shared_ptr<ComputeProgram> cp_radon_transform_fan_equiangular_parallelize_float = std::make_shared<ComputeProgram>(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/radon_transform_fan_equiangular_parallelize_float.comp"));

	std::shared_ptr<ComputeProgram> cp_int32_to_float_texture_blit		= std::make_shared<ComputeProgram>(Shader("../CTAnalyzer/Source/GLSL/Compute/Util/texture_blit_int32_to_float.comp"));
	std::shared_ptr<ComputeProgram> cp_float_to_int32_texture_blit		= std::make_shared<ComputeProgram>(Shader("../CTAnalyzer/Source/GLSL/Compute/Util/texture_blit_float_to_int32.comp"));
	std::shared_ptr<ComputeProgram> cp_float1_to_float4_texture_blit	= std::make_shared<ComputeProgram>(Shader("../CTAnalyzer/Source/GLSL/Compute/Util/texture_blit_float1_to_float4.comp"));
	std::shared_ptr<ComputeProgram> cp_float1_to_float1_texture_blit	= std::make_shared<ComputeProgram>(Shader("../CTAnalyzer/Source/GLSL/Compute/Util/texture_blit_float1_to_float1.comp"));
	std::shared_ptr<ComputeProgram> cp_texture_blit_float1_to_r16		= std::make_shared<ComputeProgram>(Shader("../CTAnalyzer/Source/GLSL/Compute/Util/texture_blit_float1_to_r16.comp"));

	std::shared_ptr<ComputeProgram> cp_int32_clip_values				= std::make_shared<ComputeProgram>(Shader("Source/GLSL/Compute/crop_overflow_pixel_values.comp"));
	std::shared_ptr<ComputeProgram> cp_clear							= std::make_shared<ComputeProgram>(Shader("../CTAnalyzer/Source/GLSL/Compute/clear.comp"));
	std::shared_ptr<ComputeProgram> cp_int_divide						= std::make_shared<ComputeProgram>(Shader("../CTAnalyzer/Source/GLSL/Compute/Util/int_texture_divide.comp"));
	std::shared_ptr<ComputeProgram> cp_create_highpass_filter			= std::make_shared<ComputeProgram>(Shader("../CTAnalyzer/Source/GLSL/Compute/create_highpass_filter.comp"));
	std::shared_ptr<ComputeProgram> cp_apply_highpass_filter			= std::make_shared<ComputeProgram>(Shader("../CTAnalyzer/Source/GLSL/Compute/apply_highpass_filter_rows.comp"));

	std::shared_ptr<ComputeProgram> cp_int16_texture_total				= std::make_shared<ComputeProgram>(Shader("../CTAnalyzer/Source/GLSL/Compute/int16_texture_total.comp"));
	std::shared_ptr<ComputeProgram> cp_int16_texture_total_2d			= std::make_shared<ComputeProgram>(Shader("../CTAnalyzer/Source/GLSL/Compute/int16_texture_total_2d.comp"));
	std::shared_ptr<ComputeProgram> cp_int_texture_subtract_2d			= std::make_shared<ComputeProgram>(Shader("../CTAnalyzer/Source/GLSL/Compute/int_texture_subtract_2d.comp"));

	std::shared_ptr<ComputeProgram> cp_fourier_space_filter_ram_lak		= std::make_shared<ComputeProgram>(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/fourier_space_filter_ram_lak.comp"));
	std::shared_ptr<ComputeProgram> cp_fourier_space_filter_shepp_logan = std::make_shared<ComputeProgram>(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/fourier_space_filter_shepp_logan.comp"));
	std::shared_ptr<ComputeProgram> cp_fourier_space_filter_cosine		= std::make_shared<ComputeProgram>(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/fourier_space_filter_cosine.comp"));

};