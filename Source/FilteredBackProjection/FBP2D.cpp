#include "FBP2D.h"
#include "FastFourierTransformCPU/FastFourierTransform.h"
#include "FastFourierTransformCPU/FastFourierTransformImage.h"

std::shared_ptr<Texture2D> FBP2D::project_forward(Texture2D& source_image, int projection_count, float total_rotation_count)
{
	std::shared_ptr<Texture2D> sinogram = std::make_shared<Texture2D>(source_image.get_size().x, projection_count, Texture2D::ColorTextureFormat::R32F, 1, 0, 0);
	project_forward(source_image, *sinogram, total_rotation_count);

	return sinogram;
}

std::shared_ptr<Texture2D> FBP2D::project_forward_equiangular(Texture2D& source_image, int projection_count, float xray_source_distance, float detector_panel_distance, float detector_width, float volume_width, float total_rotation_count)
{
	std::shared_ptr<Texture2D> sinogram = std::make_shared<Texture2D>(source_image.get_size().x, projection_count, Texture2D::ColorTextureFormat::R32F, 1, 0, 0);
	project_forward_equiangular(source_image, *sinogram, xray_source_distance, detector_panel_distance, detector_width, volume_width, total_rotation_count);

	return sinogram;
}

std::shared_ptr<Texture2D> FBP2D::project_backward(Texture2D& sinogram)
{
	glm::vec2 sinogram_size = sinogram.get_size();
	std::shared_ptr<Texture2D> slice = std::make_shared<Texture2D>(sinogram_size.x, sinogram_size.x, Texture2D::ColorTextureFormat::R32F, 1, 0, 0);

	project_backward(sinogram, *slice);

	return slice;
}

void FBP2D::parallelize_fan_beam_sinogram(Texture2D& source_sinogram, Texture2D& target_sinogram, float xray_source_distance, float detector_panel_distance, float detector_width, float volume_width, float total_rotation_count)
{
	target_sinogram.force_allocation();
	target_sinogram.clear(glm::vec4(0, 0, 0, 0), 0);

	auto& program = cp_radon_transform_fan_equiangular_parallelize_float;

	program->update_uniform_as_image("source_sinogram", source_sinogram, 0);
	program->update_uniform_as_image("target_sinogram", target_sinogram, 0);
	program->update_uniform("sinogram_resolution", target_sinogram.get_size());
	program->update_uniform("full_rotation_count", total_rotation_count);
	program->update_uniform("source_image_size", source_sinogram.get_size());
	program->update_uniform("begin_rotation_offset_radian", 0.0f);
	program->update_uniform("projection_plane_distance", detector_panel_distance);
	program->update_uniform("projection_plane_width", detector_width);
	program->update_uniform("volume_width", volume_width);
	program->update_uniform("x_ray_source_distance", xray_source_distance);

	program->update_uniform("parallel_projeciton_plane_width", volume_width * (float)std::sqrt(2));

	program->dispatch(std::ceil(source_sinogram.get_size().x / 8.0f), std::ceil(source_sinogram.get_size().y / 8.0f), 1);
}

void FBP2D::project_forward(Texture2D& source_image, Texture2D& target_image, float total_rotation_count)
{
	target_image.force_allocation();
	target_image.clear(glm::vec4(0, 0, 0, 0), 0);

	cp_radon_transform_float->update_uniform_as_image("source_image", source_image, 0);
	cp_radon_transform_float->update_uniform_as_image("target_sinogram", target_image, 0);
	cp_radon_transform_float->update_uniform("sinogram_resolution", target_image.get_size());
	cp_radon_transform_float->update_uniform("full_rotation_count", total_rotation_count);
	cp_radon_transform_float->update_uniform("source_image_size", source_image.get_size());
	cp_radon_transform_float->update_uniform("begin_rotation_offset_radian", 0.0f);
	cp_radon_transform_float->update_uniform("projection_plane_distance", 0.0f);
	cp_radon_transform_float->update_uniform("projection_plane_width", (float)source_image.get_size().x * (float)std::sqrt(2));
	cp_radon_transform_float->update_uniform("volume_width", (float)source_image.get_size().x);

	cp_radon_transform_float->dispatch(std::ceil(source_image.get_size().x / 16.0f), std::ceil(source_image.get_size().y / 16.0f), target_image.get_size().y);

	//_int_texture_div(target_image, target_image.get_size().x);
}

void FBP2D::project_forward_equiangular(Texture2D& source_image, Texture2D& target_image, float xray_source_distance, float detector_panel_distance, float detector_width, float volume_width, float total_rotation_count)
{
	target_image.force_allocation();
	target_image.clear(glm::vec4(0, 0, 0, 0), 0);

	cp_radon_transform_fan_equiangular_float->update_uniform_as_image("source_image", source_image, 0);
	cp_radon_transform_fan_equiangular_float->update_uniform_as_image("target_sinogram", target_image, 0);
	cp_radon_transform_fan_equiangular_float->update_uniform("sinogram_resolution", target_image.get_size());
	cp_radon_transform_fan_equiangular_float->update_uniform("full_rotation_count", total_rotation_count);
	cp_radon_transform_fan_equiangular_float->update_uniform("source_image_size", source_image.get_size());
	cp_radon_transform_fan_equiangular_float->update_uniform("begin_rotation_offset_radian", 0.0f);
	cp_radon_transform_fan_equiangular_float->update_uniform("x_ray_source_distance", xray_source_distance);
	cp_radon_transform_fan_equiangular_float->update_uniform("projection_plane_distance", detector_panel_distance);
	cp_radon_transform_fan_equiangular_float->update_uniform("projection_plane_width", detector_width);
	cp_radon_transform_fan_equiangular_float->update_uniform("volume_width", volume_width);


	cp_radon_transform_fan_equiangular_float->dispatch(std::ceil(source_image.get_size().x / 16.0f), std::ceil(source_image.get_size().y / 16.0f), target_image.get_size().y);

	//_int_texture_div(target_image, target_image.get_size().x);
}

void FBP2D::project_backward(Texture2D& sinogram, Texture2D& out_texture)
{
	out_texture.force_allocation();

	cp_clear->update_uniform_as_image("target", out_texture, 0);
	cp_clear->dispatch(std::ceil(out_texture.get_size().x / 8.0f), std::ceil(out_texture.get_size().y / 8.0f), 1);

	cp_radon_transform_inverse_float->update_uniform_as_image("sinogram", sinogram, 0);
	cp_radon_transform_inverse_float->update_uniform_as_image("target_slice", out_texture, 0);
	cp_radon_transform_inverse_float->update_uniform("sinogram_resolution", (int)sinogram.get_size().x, (int)sinogram.get_size().y);
	cp_radon_transform_inverse_float->update_uniform("slice_resolution", (int)out_texture.get_size().x, (int)out_texture.get_size().y);
	cp_radon_transform_inverse_float->update_uniform("full_rotation_count", 1.0f);
	cp_radon_transform_inverse_float->update_uniform("begin_rotation_offset_radian", 0.006f);
	cp_radon_transform_inverse_float->update_uniform("projection_plane_distance", 0.0f);
	cp_radon_transform_inverse_float->update_uniform("projection_plane_width", (float)out_texture.get_size().x * (float)std::sqrt(2));
	cp_radon_transform_inverse_float->update_uniform("volume_width", (float)out_texture.get_size().x);

	cp_radon_transform_inverse_float->dispatch(std::ceil(out_texture.get_size().x / 16.0f), std::ceil(out_texture.get_size().y / 16.0f), std::ceil(sinogram.get_size().y / 1.0f));

	//unsigned int sinogram_total_avarage = _total_value_of_texture(sinogram) / sinogram.get_size().y;
	//std::cout << "total is: " << sinogram_total_avarage << std::endl;
	//_subtract_from_texture(out_texture, sinogram_total_avarage);

	cp_int_divide->update_uniform_as_image("int_texture", out_texture, 0);
	cp_int_divide->update_uniform("denominator", (float)sinogram.get_size().y);
	cp_int_divide->dispatch(std::ceil(out_texture.get_size().x / 8.0f), std::ceil(out_texture.get_size().y / 8.0f), 1);
}

void FBP2D::apply_filter(Texture2D& sinogram, FilterType filter_type)
{
	//std::shared_ptr<Image> sinogram_image = source_sinogram.get_image(Texture2D::ColorFormat::RED, Texture2D::Type::UNSIGNED_INT, 0);
	//FFT::image_fft_1d_radix2(*sinogram_image);
	//target_sinogram.load_data(*sinogram_image, Texture2D::ColorFormat::RED, Texture2D::Type::UNSIGNED_INT, 0);

	std::shared_ptr<ComputeProgram> program_to_use = 
		filter_type == FilterType::RAM_LAK		? cp_fourier_space_filter_ram_lak		:
		filter_type == FilterType::SHEPP_LOGAN	? cp_fourier_space_filter_shepp_logan	:
		filter_type == FilterType::COSINE		? cp_fourier_space_filter_cosine : nullptr;
	
	if (program_to_use == nullptr) {
		ASSERT(false);
	}

	program_to_use->update_uniform_as_image("fourier_texture", sinogram, 0);
	program_to_use->update_uniform("texture_resolution", sinogram.get_size());

	program_to_use->dispatch(std::ceil(sinogram.get_size().x / 8.0f), std::ceil(sinogram.get_size().y / 8.0f), 1);

}

void FBP2D::_texture_blit_int32_to_float(Texture2D& source_int_texture, Texture2D& target_float_texture)
{
	cp_int32_to_float_texture_blit->update_uniform_as_image("int_texture", source_int_texture, 0);
	cp_int32_to_float_texture_blit->update_uniform_as_image("target_texture", target_float_texture, 0);
	cp_int32_to_float_texture_blit->dispatch(std::ceil(target_float_texture.get_size().x / 8.0f), std::ceil(target_float_texture.get_size().y / 8.0f), 1);
}

void FBP2D::_texture_blit_float_to_int32(Texture2D& source_float_texture, Texture2D& target_int_texture)
{
	cp_float_to_int32_texture_blit->update_uniform_as_image("int_texture", source_float_texture, 0);
	cp_float_to_int32_texture_blit->update_uniform_as_image("target_texture", target_int_texture, 0);
	cp_float_to_int32_texture_blit->dispatch(std::ceil(target_int_texture.get_size().x / 8.0f), std::ceil(target_int_texture.get_size().y / 8.0f), 1);
}

void FBP2D::_texture_blit_float1_to_float4(Texture2D& source_float_texture, Texture2D& target_float_texture)
{
	cp_float1_to_float4_texture_blit->update_uniform_as_image("source_float_texture", source_float_texture, 0);
	cp_float1_to_float4_texture_blit->update_uniform_as_image("target_float_texture", target_float_texture, 0);
	cp_float1_to_float4_texture_blit->dispatch(std::ceil(target_float_texture.get_size().x / 8.0f), std::ceil(target_float_texture.get_size().y / 8.0f), 1);

}

void FBP2D::_texture_blit_float1_to_float1(Texture2D& source_float_texture, Texture2D& target_float_texture)
{
	cp_float1_to_float1_texture_blit->update_uniform_as_image("source_float_texture", source_float_texture, 0);
	cp_float1_to_float1_texture_blit->update_uniform_as_image("target_float_texture", target_float_texture, 0);
	cp_float1_to_float1_texture_blit->dispatch(std::ceil(target_float_texture.get_size().x / 8.0f), std::ceil(target_float_texture.get_size().y / 8.0f), 1);
}

void FBP2D::_texture_blit_float1_to_r16(Texture2D& source_float_texture, Texture2D& target_float_texture)
{
	cp_texture_blit_float1_to_r16->update_uniform_as_image("source_float_texture", source_float_texture, 0);
	cp_texture_blit_float1_to_r16->update_uniform_as_image("target_float_texture", target_float_texture, 0);
	cp_texture_blit_float1_to_r16->dispatch(std::ceil(target_float_texture.get_size().x / 8.0f), std::ceil(target_float_texture.get_size().y / 8.0f), 1);
}

void FBP2D::_int_texture_clip_overflow(Texture2D& int_texture)
{
	cp_int32_clip_values->update_uniform_as_image("int_texture", int_texture, 0);
	cp_int32_clip_values->update_uniform("image_resolution", int_texture.get_size());
	cp_int32_clip_values->dispatch(std::ceil(int_texture.get_size().x / 8.0f), std::ceil(int_texture.get_size().y / 8.0f), 1);
}

void FBP2D::_int_texture_div(Texture2D& int_texture, float denominator)
{
	cp_int_divide->update_uniform_as_image("int_texture", int_texture, 0);
	cp_int_divide->update_uniform("texture_resolution", int_texture.get_size());
	cp_int_divide->update_uniform("denominator", denominator);
	cp_int_divide->dispatch(std::ceil(int_texture.get_size().x / 8.0f), std::ceil(int_texture.get_size().y / 8.0f), 1);
}
