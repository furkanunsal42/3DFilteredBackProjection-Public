#include "FBP3D.h"
#include <fstream>
#include <iostream>
#include <cstdint>
#include <filesystem>
#include <mutex>
#include <chrono>

#include "ComputeProgram.h"

FBP3D::FBP3D()
{

}

void FBP3D::generate_blank_volume(int x_size, int y_size, int z_size)
{
	if (volume != nullptr) {
		if (volume->get_size() == glm::ivec3(x_size, y_size, z_size)) {
			volume->force_allocation();
			volume->clear(0.0f, 0);
			return;
		}
	}

	volume = std::make_shared<Texture3D>(x_size, y_size, z_size, _volume_internal_format, 1, 0);
	volume->is_bindless = false;
}

void FBP3D::generate_blank_projections(int x_size, int y_size, int projection_amount)
{
	if (projections != nullptr) {
		if (projections->get_size() == glm::ivec3(x_size, y_size, projection_amount)) {
			projections->force_allocation();
			projections->clear(0.0f, 0);
			return;
		}
	}

	projections = std::make_shared<Texture2DArray>(x_size, y_size, projection_amount, _projection_internal_format, 1, 0, 0);
	projections->is_bindless = false;
}

void FBP3D::generate_blank_sinograms(int x_size, int y_size, int projection_amount)
{
	if (sinograms != nullptr) {
		if (sinograms->get_size() == glm::ivec3(x_size, projection_amount, y_size)) {
			sinograms->force_allocation();
			sinograms->clear(0.0f, 0);
			return;
		}
	}

	sinograms = std::make_shared<Texture3D>(x_size, projection_amount, y_size, _projection_internal_format, 1, 0);
	sinograms->is_bindless = false;
}

void FBP3D::generate_shepplogan(int x_size, int y_size, int z_size)
{
	generate_blank_volume(x_size, y_size, z_size);
	
	render_ellipsoid(glm::vec3(0	 , 0.0f	  , 0.0f   ), glm::vec3(0.69  , 0.9  , 0.92 ), glm::vec3(0, 0  , 0) / 180.0f * 3.14f, (2.0  ) / 3.0f);
	render_ellipsoid(glm::vec3(0	 , 0.0f	  , 0.0f   ), glm::vec3(0.6624, 0.88 , 0.874), glm::vec3(0, 0  , 0) / 180.0f * 3.14f, (-0.98) / 3.0f);
	render_ellipsoid(glm::vec3(-0.22f, -0.25f , 0.0f   ), glm::vec3(0.41  , 0.21 , 0.16 ), glm::vec3(0, 108, 0) / 180.0f * 3.14f, (-0.02) * 75 / 3.0f);
	render_ellipsoid(glm::vec3(0.22f , -0.25f , 0.0f   ), glm::vec3(0.31  , 0.22 , 0.11 ), glm::vec3(0, 72 , 0) / 180.0f * 3.14f, (-0.02) * 75 / 3.0f);
	render_ellipsoid(glm::vec3(0	 , -0.25f , 0.4f   ), glm::vec3(0.23  , 0.46 , 0.23 ), glm::vec3(0, 0  , 0) / 180.0f * 3.14f, (0.02 ) * 75 / 3.0f);
	render_ellipsoid(glm::vec3(0	 , -0.25f , 0.1f   ), glm::vec3(0.046 , 0.046, 0.046), glm::vec3(0, 0  , 0) / 180.0f * 3.14f, (0.02 ) * 75 / 3.0f);
	render_ellipsoid(glm::vec3(-0.8f , -0.25f , -0.65f ), glm::vec3(0.046 , 0.02 , 0.023), glm::vec3(0, 0  , 0) / 180.0f * 3.14f, (0.01 ) * 75 / 3.0f);
	render_ellipsoid(glm::vec3(0.06f , -0.25f , -0.065f), glm::vec3(0.046 , 0.02 , 0.023), glm::vec3(0, 90 , 0) / 180.0f * 3.14f, (0.01 ) * 75 / 3.0f);
	render_ellipsoid(glm::vec3(0.06f , 0.625f , -0.105f), glm::vec3(0.56  , 0.1  , 0.04 ), glm::vec3(0, 90 , 0) / 180.0f * 3.14f, (0.02 ) * 75 / 3.0f);
	render_ellipsoid(glm::vec3(0	 , -0.625f, 0.1f   ), glm::vec3(0.056 , 0.1  , 0.056), glm::vec3(0, 0  , 0) / 180.0f * 3.14f, (-0.02) * 75 / 3.0f);
}

void FBP3D::read_projections(const std::string& upper_diretory, int file_width, int file_height, int file_channel_count, int byte_per_channel, int target_projection_width, int target_projection_height, int projection_count)
{
	generate_blank_projections(target_projection_width, target_projection_height, projection_count);

	std::shared_ptr<Texture2D> temp_projection = std::make_shared<Texture2D>(target_projection_width, target_projection_height, _projection_internal_format, 1, 0, 0);
	temp_projection->is_bindless = false;

	std::cout << "[FBP3D Info] loading from \"" << upper_diretory << "\"" << std::endl;
	auto begin_time = std::chrono::system_clock::now();
	const bool print_info = false;

	std::vector<std::shared_ptr<Image>> ct_projections(projection_count);
	std::vector<std::shared_ptr<std::thread>> loading_threads(projection_count);

	if (!std::filesystem::exists(upper_diretory)) {
		std::cout << "[CTAnalyzer Error] FBP3D tried to load_projections() but \"" << upper_diretory << "\" doesn't exist" << std::endl;
		ASSERT(false);
	}
	else {

		std::vector<std::string> files_in_directory(projection_count);
		for (const auto& path : std::filesystem::directory_iterator(upper_diretory)) {
			std::string& filepath = path.path().string();
			if (filepath.size() >= 8) {
				int slice_index = std::stoi(filepath.substr(filepath.size() - 8, 4));
				if (slice_index >= projection_count) continue;
				if (slice_index < 0) continue;
				files_in_directory[slice_index] = filepath;
			}
		}

		const int batch_size = 128;
		for (int batch_index = 0; batch_size * batch_index < projection_count; batch_index++) {
			#pragma omp parallel for num_threads(16) shared(ct_projections)
			for (int i = batch_index * batch_size; i < (int)std::min((batch_index + 1) * batch_size, projection_count); i++) {
				ct_projections[i] = std::make_shared<Image>(files_in_directory[i], file_width, file_height, 1, file_channel_count, byte_per_channel, true);
				ct_projections[i]->resize(target_projection_width, target_projection_height);
			}

			for (int i = batch_index * batch_size; i < (int)std::min((batch_index + 1) * batch_size, projection_count); i++) {
				temp_projection->load_data(*ct_projections[i], Texture2D::ColorFormat::RED, Texture2D::Type::UNSIGNED_SHORT, 0);
				store_projection(*temp_projection, i);
				ct_projections[i] = nullptr;
			}
		}
	}

	auto end_time = std::chrono::system_clock::now();
	std::cout << "[FBP3D Info] load finished in " << std::chrono::duration<double>(end_time - begin_time).count() << " seconds" << std::endl;

}

void FBP3D::read_projections_as_sinograms(const std::string& upper_diretory, int file_width, int file_height, int file_channel_count, int byte_per_channel, int target_projection_width, int target_projection_height, int projection_count)
{
	generate_blank_sinograms(target_projection_width, target_projection_height, projection_count);

	std::shared_ptr<Texture2D> temp_projection = std::make_shared<Texture2D>(target_projection_width, target_projection_height, _projection_internal_format, 1, 0, 0);
	temp_projection->is_bindless = false;

	std::cout << "[FBP3D Info] loading from \"" << upper_diretory << "\"" << std::endl;
	auto begin_time = std::chrono::system_clock::now();
	const bool print_info = false;

	std::vector<std::shared_ptr<Image>> ct_projections(projection_count);
	std::vector<std::shared_ptr<std::thread>> loading_threads(projection_count);

	if (!std::filesystem::exists(upper_diretory)) {
		std::cout << "[CTAnalyzer Error] FBP3D tried to load_sinograms() but \"" << upper_diretory << "\" doesn't exist" << std::endl;
		ASSERT(false);
	}
	else {

		std::vector<std::string> files_in_directory(projection_count);
		for (const auto& path : std::filesystem::directory_iterator(upper_diretory)) {
			std::string& filepath = path.path().string();
			if (filepath.size() >= 8) {
				int slice_index = std::stoi(filepath.substr(filepath.size() - 8, 4));
				if (slice_index >= projection_count) continue;
				if (slice_index < 0) continue;
				files_in_directory[slice_index] = filepath;
			}
		}

		const int batch_size = 128;
		for (int batch_index = 0; batch_size * batch_index < projection_count; batch_index++) {
			#pragma omp parallel for num_threads(16) shared(ct_projections)
			for (int i = batch_index * batch_size; i < (int)std::min((batch_index + 1) * batch_size, projection_count); i++) {
				ct_projections[i] = std::make_shared<Image>(files_in_directory[i], file_width, file_height, 1, file_channel_count, byte_per_channel, true);
				ct_projections[i]->resize(target_projection_width, target_projection_height);
			}

			for (int i = batch_index * batch_size; i < (int)std::min((batch_index + 1) * batch_size, projection_count); i++) {
				temp_projection->load_data(*ct_projections[i], Texture2D::ColorFormat::RED, Texture2D::Type::UNSIGNED_SHORT, 0);
				store_projection_to_sinograms(*temp_projection, i);
				ct_projections[i] = nullptr;
			}
		}
	}

	auto end_time = std::chrono::system_clock::now();
	std::cout << "[FBP3D Info] load finished in " << std::chrono::duration<double>(end_time - begin_time).count() << " seconds" << std::endl;

}

void FBP3D::write_volume(std::string upper_directory)
{
	std::shared_ptr<Texture2D> slice = std::make_shared<Texture2D>(volume->get_size().x, volume->get_size().y, _volume_internal_format, 1, 0, 0);
	slice->is_bindless = false;

	if (upper_directory.at(upper_directory.size() - 1) != '/' || upper_directory.at(upper_directory.size() - 1) != '\\') upper_directory.insert(upper_directory.end(), '/');
	
	for (int i = 0; i < volume->get_size().y; i++) {
		
		load_volume_slice_y(i, *slice);

		std::string slice_index = std::to_string(i);
		slice_index.insert(slice_index.begin(), _slice_index_number_length - slice_index.size(), '0');
		std::shared_ptr<Image> slice_image = slice->get_image(Texture2D::ColorFormat::RED, Texture2D::Type::UNSIGNED_SHORT, 0);
		slice_image->save_to_disc(upper_directory + "CTAnalyzer_Slice#" + slice_index + ".raw");

	}
}

void FBP3D::free_volume()
{
	volume = nullptr;
}

void FBP3D::free_projections()
{
	projections = nullptr;
}

void FBP3D::free_sinograms()
{
	sinograms = nullptr;
}

void FBP3D::project_forward_parallel_to_projections(float full_rotation_count, int projection_count, float rotation_offset_radian)
{
	if (volume == nullptr) {
		std::cout << "[CTAnalyzer Error] FBP3D tried to project_forward_parallel() but FBP3D.volume was nullptr" << std::endl;
		ASSERT(false);
	}

	generate_blank_projections(volume->get_size().x, volume->get_size().y, projection_count);

	ComputeProgram& kernel =  *cp_ct3d_radon_transform_float_parallel;

	const float PI = 3.14159265359f;
	float total_rotation = 2 * PI * full_rotation_count;
	for (int invocation_id = 0; invocation_id < projection_count; invocation_id++) {
		float invocation_angle_radian = total_rotation / projection_count * invocation_id;
		invocation_angle_radian += rotation_offset_radian;

		kernel.update_uniform_as_image("source_volume", *volume, 0);
		kernel.update_uniform_as_image("target_projections", *projections, 0);
		kernel.update_uniform("volume_resolution", volume->get_size());
		kernel.update_uniform("projection_resolution", projections->get_size().x, projections->get_size().y);
		kernel.update_uniform("projection_plane_width", (float)(volume->get_size().x));
		kernel.update_uniform("volume_width", (float)volume->get_size().x);
		kernel.update_uniform("projection_vector", glm::vec2(glm::cos(invocation_angle_radian), glm::sin(invocation_angle_radian)));
		kernel.update_uniform("projection_index", invocation_id);

		kernel.dispatch(std::ceil(projections->get_size().x / 8.0f), std::ceil(projections->get_size().y / 8.0f), std::ceil(volume->get_size().x/ 1));
	}
}

void FBP3D::project_forward_parallel_to_sinogram(float full_rotation_count, int projection_count, float rotation_offset_radian)
{
	if (volume == nullptr) {
		std::cout << "[CTAnalyzer Error] FBP3D tried to project_forward_parallel_to_sinogram() but FBP3D.volume was nullptr" << std::endl;
		ASSERT(false);
	}

	generate_blank_sinograms(volume->get_size().x, volume->get_size().y, projection_count);

	ComputeProgram& kernel = *cp_ct3d_radon_transform_float_parallel_sinogram;

	int height = volume->get_size().y;
	for (int invocation_id = 0; invocation_id < height; invocation_id++) {
		kernel.update_uniform("source_volume", *volume);
		kernel.update_uniform_as_image("target_sinograms", *sinograms, 0);
		kernel.update_uniform("volume_resolution", volume->get_size());
		kernel.update_uniform("sinogram_resolution", sinograms->get_size().x, sinograms->get_size().y);
		kernel.update_uniform("projection_plane_width", (float)(volume->get_size().x));
		kernel.update_uniform("volume_width", (float)volume->get_size().x);
		kernel.update_uniform("full_rotation_count", full_rotation_count);
		kernel.update_uniform("rotation_offset_radian", rotation_offset_radian);
		kernel.update_uniform("slice_y_index", invocation_id);

		kernel.dispatch(std::ceil(sinograms->get_size().x / 8.0f), std::ceil(sinograms->get_size().y / 8.0f), std::ceil(volume->get_size().x / 1));
	}
}

void FBP3D::project_forward_cone_to_projections(float xray_source_distance, float detector_panel_distance, float detector_width, float volume_width, float full_rotation_count, int projection_count, float rotation_offset_radian)
{
	if (volume == nullptr) {
		std::cout << "[CTAnalyzer Error] FBP3D tried to project_forward_cone_to_projections() but FBP3D.volume was nullptr" << std::endl;
		ASSERT(false);
	}

	generate_blank_projections(volume->get_size().x, volume->get_size().y, projection_count);

	ComputeProgram& kernel = *cp_ct3d_radon_transform_float_cone_equidistance;
	
	const float PI = 3.14159265359f;
	float total_rotation = 2 * PI * full_rotation_count;
	for (int invocation_id = 0; invocation_id < projection_count; invocation_id++) {
		float invocation_angle_radian = total_rotation / projection_count * invocation_id;
		invocation_angle_radian += rotation_offset_radian;

		kernel.update_uniform("source_volume", *volume);
		kernel.update_uniform_as_image("target_projections", *projections, 0);
		kernel.update_uniform("volume_resolution", volume->get_size());
		kernel.update_uniform("projection_resolution", projections->get_size().x, projections->get_size().y);
		kernel.update_uniform("projection_plane_width", detector_width);
		kernel.update_uniform("volume_width", volume_width);
		kernel.update_uniform("center_to_detector_distance", detector_panel_distance);
		kernel.update_uniform("center_to_source_distance", xray_source_distance);
		kernel.update_uniform("projection_vector", glm::vec2(glm::cos(invocation_angle_radian), sin(invocation_angle_radian)));
		kernel.update_uniform("projection_index", invocation_id);

		kernel.dispatch(std::ceil(projections->get_size().x / 8.0f), std::ceil(projections->get_size().y / 8.0f), std::ceil(volume->get_size().x / 1));
	}
}

void FBP3D::project_forward_cone_to_sinogram(float xray_source_distance, float detector_panel_distance, float detector_width, float volume_width, float full_rotation_count, int projection_count, float rotation_offset_radian)
{
	if (volume == nullptr) {
		std::cout << "[CTAnalyzer Error] FBP3D tried to project_forward_cone_to_sinogram() but FBP3D.volume was nullptr" << std::endl;
		ASSERT(false);
	}

	generate_blank_sinograms(volume->get_size().x, volume->get_size().y, projection_count);

	ComputeProgram& kernel = *cp_ct3d_radon_transform_float_cone_equidistance_sinogram;

	int height = sinograms->get_size().z;
	for (int invocation_id = 0; invocation_id < height ; invocation_id++) {

		kernel.update_uniform("source_volume", *volume);
		kernel.update_uniform_as_image("target_sinograms", *sinograms, 0);
		kernel.update_uniform("volume_resolution", volume->get_size());
		kernel.update_uniform("sinogram_resolution", sinograms->get_size().x, sinograms->get_size().y);
		kernel.update_uniform("sinogram_count", sinograms->get_size().z);
		kernel.update_uniform("projection_plane_width", detector_width);
		kernel.update_uniform("volume_width", volume_width);
		kernel.update_uniform("center_to_detector_distance", detector_panel_distance);
		kernel.update_uniform("center_to_source_distance", xray_source_distance);
		kernel.update_uniform("full_rotation_count", full_rotation_count);
		kernel.update_uniform("rotation_offset_radian", rotation_offset_radian);
		kernel.update_uniform("slice_y_index", invocation_id);

		kernel.dispatch(std::ceil(sinograms->get_size().x / 8.0f), std::ceil(sinograms->get_size().y / 8.0f), std::ceil(volume->get_size().x / 1));
	}
}


void FBP3D::log_normalize_projections(float base_value)
{
	std::shared_ptr<Texture2D> projection = std::make_shared<Texture2D>(projections->get_size().x, projections->get_size().y, _projection_internal_format, 1, 0, 0);
	projection->is_bindless = false;

	for (int projection_index = 0; projection_index < projections->get_size().z; projection_index++) {
		load_projection(projection_index, *projection);

		cp_ct3d_log_normalize->update_uniform_as_image("source_texture", *projection, 0);
		cp_ct3d_log_normalize->update_uniform("source_texture_resolution", projection->get_size());
		cp_ct3d_log_normalize->update_uniform("base_value", base_value);

		cp_ct3d_log_normalize->dispatch(std::ceil(projection->get_size().x / 8.0f), std::ceil(projection->get_size().y / 8.0f), 1);
		
		store_projection(*projection, projection_index);
	}
}

void FBP3D::log_normalize_sinograms(float base_value)
{
	std::shared_ptr<Texture2D> sinogram = std::make_shared<Texture2D>(projections->get_size().x, projections->get_size().y, _projection_internal_format, 1, 0, 0);
	sinogram->is_bindless = false;

	for (int projection_index = 0; projection_index < projections->get_size().z; projection_index++) {
		load_sinogram(projection_index, *sinogram);

		cp_ct3d_log_normalize->update_uniform_as_image("source_texture", *sinogram, 0);
		cp_ct3d_log_normalize->update_uniform("source_texture_resolution", sinogram->get_size());
		cp_ct3d_log_normalize->update_uniform("base_value", base_value);

		cp_ct3d_log_normalize->dispatch(std::ceil(sinogram->get_size().x / 8.0f), std::ceil(sinogram->get_size().y / 8.0f), 1);

		store_sinogram(*sinogram, projection_index);
	}
}

void FBP3D::apply_fdk_weights_to_projections(float xray_source_distance, float detector_panel_distance, float detector_width)
{
	std::shared_ptr<Texture2D> projection = std::make_shared<Texture2D>(projections->get_size().x, projections->get_size().y, _projection_internal_format, 1, 0, 0);
	projection->is_bindless = false;

	for (int projection_index = 0; projection_index < projections->get_size().z; projection_index++) {
		load_projection(projection_index, *projection);

		cp_ct3d_fdk_weight->update_uniform_as_image("source_texture", *projection, 0);
		cp_ct3d_fdk_weight->update_uniform("source_texture_resolution", projection->get_size());
		cp_ct3d_fdk_weight->update_uniform("projection_plane_width", detector_width);
		cp_ct3d_fdk_weight->update_uniform("center_to_detector_distance", detector_panel_distance);
		cp_ct3d_fdk_weight->update_uniform("center_to_source_distance", xray_source_distance);

		cp_ct3d_fdk_weight->dispatch(std::ceil(projection->get_size().x / 8.0f), std::ceil(projection->get_size().y / 8.0f), 1);

		store_projection(*projection, projection_index);
	}
}

void FBP3D::apply_filter_to_projections(FBP2D::FilterType filter_type)
{
	std::shared_ptr<Texture2D> projection = std::make_shared<Texture2D>(projections->get_size().x, projections->get_size().y, _projection_internal_format, 1, 0, 0);
	projection->is_bindless = false;
	std::shared_ptr<Texture2D> projection_complex = std::make_shared<Texture2D>(projections->get_size().x, projections->get_size().y, _projection_complex_internal_format, 1, 0, 0);
	projection_complex->is_bindless = false;

	std::shared_ptr<Texture2D> projection_complex_padded = fft_solver->created_padded_complex_texture(*projection_complex);
	std::cout << projection_complex_padded->get_size().x << " " << projection_complex_padded->get_size().y << std::endl;

	for (int projection_index = 0; projection_index < projections->get_size().z; projection_index++) {
		load_projection(projection_index, *projection);

		fft_solver->blit_texture_real_to_complex(*projection, *projection_complex);
		fft_solver->blit_texture_complex_to_complex(*projection_complex, *projection_complex_padded);

		fft_solver->fft(*projection_complex, *projection_complex_padded);
		fft_solver->fft_shift(*projection_complex_padded);
		fbp_solver->apply_filter(*projection_complex_padded, filter_type);
		fft_solver->inverse_fft_shift(*projection_complex_padded);
		fft_solver->inverse_fft(*projection_complex_padded, *projection_complex);
	
		fft_solver->blit_texture_complex_to_real_r(*projection_complex, *projection);
		store_projection(*projection, projection_index);
	}
}


void FBP3D::apply_filter_to_sinograms(FBP2D::FilterType filter_type)
{
	std::shared_ptr<Texture2D> sinogram = std::make_shared<Texture2D>(sinograms->get_size().x, sinograms->get_size().y, _projection_internal_format, 1, 0, 0);
	sinogram->is_bindless = false;
	std::shared_ptr<Texture2D> sinogram_complex = std::make_shared<Texture2D>(sinograms->get_size().x, sinograms->get_size().y, _projection_complex_internal_format, 1, 0, 0);
	sinogram->is_bindless = false;

	std::shared_ptr<Texture2D> sinogram_complex_padded = fft_solver->created_padded_complex_texture(*sinogram_complex);
	std::cout << sinogram_complex_padded->get_size().x << " " << sinogram_complex_padded->get_size().y << std::endl;

	for (int sinogram_index = 0; sinogram_index < sinograms->get_size().z; sinogram_index++) {
		load_sinogram(sinogram_index, *sinogram);

		fft_solver->blit_texture_real_to_complex(*sinogram, *sinogram_complex);
		fft_solver->blit_texture_complex_to_complex(*sinogram_complex, *sinogram_complex_padded);

		fft_solver->fft(*sinogram_complex, *sinogram_complex_padded);
		fft_solver->fft_shift(*sinogram_complex_padded);
		fbp_solver->apply_filter(*sinogram_complex_padded, filter_type);
		fft_solver->inverse_fft_shift(*sinogram_complex_padded);
		fft_solver->inverse_fft(*sinogram_complex_padded, *sinogram_complex);

		fft_solver->blit_texture_complex_to_real_r(*sinogram_complex, *sinogram);
		store_sinogram(*sinogram, sinogram_index);
	}
}

void FBP3D::clip_negatives_of_volume()
{
	std::shared_ptr<Texture2D> slice = std::make_shared<Texture2D>(volume->get_size().x, volume->get_size().y, _volume_internal_format, 1, 0, 0);
	slice->is_bindless = false;
	for (int slice_index = 0; slice_index < volume->get_size().z; slice_index++) {
		load_volume_slice_y(slice_index, *slice);

		cp_ct3d_clip_negative->update_uniform_as_image("source_texture", *slice, 0);
		cp_ct3d_clip_negative->update_uniform("source_texture_resolution", slice->get_size());

		cp_ct3d_clip_negative->dispatch(std::ceil(slice->get_size().x / 8.0f), std::ceil(slice->get_size().y / 8.0f), 1);
		
		store_volume_slice_y(*slice, slice_index);
	}
}

void FBP3D::project_backward_parallel_from_projections(float full_rotation_count, int volume_resolution_horizontal, int volume_resolution_vertical, float rotation_offset_radian)
{
	if (projections == nullptr) {
		std::cout << "[CTAnalyzer Error] FBP3D tried to project_backward_parallel_from_projections() but FBP3D.projections was nullptr" << std::endl;
		ASSERT(false);
	}

	generate_blank_volume(volume_resolution_horizontal, volume_resolution_vertical, volume_resolution_horizontal);

	ComputeProgram& kernel = *cp_ct3d_radon_transform_inverse_float_parallel;

	const float PI = 3.14159265359f;
	float total_rotation = 2 * PI * full_rotation_count;
	for (int invocation_id = 0; invocation_id < projections->get_size().z; invocation_id++) {
		float invocation_angle_radian = total_rotation / projections->get_size().z * invocation_id;
		invocation_angle_radian += rotation_offset_radian;

		kernel.update_uniform_as_image("source_projections", *projections, 0);
		kernel.update_uniform_as_image("target_volume", *volume, 0);
		kernel.update_uniform("volume_resolution", volume->get_size());
		kernel.update_uniform("projection_resolution", projections->get_size().x, projections->get_size().y);
		kernel.update_uniform("projection_plane_width", (float)(volume->get_size().x));
		kernel.update_uniform("volume_width", (float)volume->get_size().x);
		kernel.update_uniform("projection_vector", glm::vec2(glm::cos(invocation_angle_radian), sin(invocation_angle_radian)));
		kernel.update_uniform("projection_index", invocation_id);

		kernel.dispatch(std::ceil(projections->get_size().x / 8.0f), std::ceil(projections->get_size().y / 8.0f), std::ceil(volume->get_size().x / 1));
	}

	clip_negatives_of_volume();
}

void FBP3D::project_backward_parallel_from_sinograms(float full_rotation_count, int volume_resolution_horizontal, int volume_resolution_vertical, float rotation_offset_radian)
{
	if (sinograms == nullptr) {
		std::cout << "[CTAnalyzer Error] FBP3D tried to project_backward_parallel_from_sinograms() but FBP3D.sinograms was nullptr" << std::endl;
		ASSERT(false);
	}

	generate_blank_volume(volume_resolution_horizontal, volume_resolution_vertical, volume_resolution_horizontal);

	ComputeProgram& kernel = *cp_ct3d_radon_transform_inverse_float_parallel_sinogram;

	int height = volume->get_size().y;
	for (int invocation_id = 0; invocation_id < height; invocation_id++) {
		
		kernel.update_uniform("source_sinograms", *sinograms);
		kernel.update_uniform_as_image("target_volume", *volume, 0);
		kernel.update_uniform("volume_resolution", volume->get_size());
		kernel.update_uniform("sinogram_resolution", sinograms->get_size().x, sinograms->get_size().y);
		kernel.update_uniform("sinogram_count", sinograms->get_size().z);
		kernel.update_uniform("projection_plane_width", (float)(volume->get_size().x));
		kernel.update_uniform("volume_width", (float)volume->get_size().x);

		kernel.update_uniform("full_rotation_count", full_rotation_count);
		kernel.update_uniform("rotation_offset_radian", rotation_offset_radian);
		kernel.update_uniform("slice_y_index", invocation_id);

		kernel.dispatch(std::ceil(volume->get_size().x / 8.0f), std::ceil(volume->get_size().y / 8.0f), std::ceil(sinograms->get_size().y / 1));
	}

	clip_negatives_of_volume();
}

void FBP3D::project_backward_cone_fdk_from_projections(float xray_source_distance, float detector_panel_distance, float detector_width, float volume_width, float full_rotation_count, int volume_resolution_horizontal, int volume_resolution_vertical, float rotation_offset_radian)
{
	if (projections == nullptr) {
		std::cout << "[CTAnalyzer Error] FBP3D tried to project_backward_cone_fdk() but FBP3D.projections was nullptr" << std::endl;
		ASSERT(false);
	}

	generate_blank_volume(volume_resolution_horizontal, volume_resolution_vertical, volume_resolution_horizontal);

	//ComputeProgram& kernel = *cp_ct3d_radon_transform_inverse_float_cone_equidistance;
	ComputeProgram& kernel = *cp_ct3d_radon_transform_inverse_float_cone_equidistance_fast;

	const float PI = 3.14159265359f;
	float total_rotation = 2 * PI * full_rotation_count;
	for (int invocation_id = 0; invocation_id < projections->get_size().z; invocation_id++) {
		float invocation_angle_radian = total_rotation / projections->get_size().z * invocation_id;
		invocation_angle_radian += rotation_offset_radian;

		//kernel.update_uniform_as_image("source_projections", *projections, 0);
		kernel.update_uniform("source_projections", *projections);
		kernel.update_uniform_as_image("target_volume", *volume, 0);
		kernel.update_uniform("volume_resolution", volume->get_size());
		kernel.update_uniform("projection_resolution", projections->get_size().x, projections->get_size().y);
		kernel.update_uniform("projection_count", projections->get_size().z);
		kernel.update_uniform("projection_plane_width", detector_width);
		kernel.update_uniform("volume_width", volume_width);
		kernel.update_uniform("center_to_detector_distance", detector_panel_distance);
		kernel.update_uniform("center_to_source_distance", xray_source_distance);
		kernel.update_uniform("projection_vector", glm::vec2(glm::cos(invocation_angle_radian), sin(invocation_angle_radian)));
		kernel.update_uniform("projection_index", invocation_id);

		kernel.dispatch(std::ceil(projections->get_size().x / 8.0f), std::ceil(projections->get_size().y / 8.0f), std::ceil(volume->get_size().x / 1));
	}

	clip_negatives_of_volume();
}

void FBP3D::project_backward_cone_fdk_from_sinograms(float xray_source_distance, float detector_panel_distance, float detector_width, float volume_width, float full_rotation_count, int volume_resolution_horizontal, int volume_resolution_vertical, float rotation_offset_radian)
{
	if (sinograms == nullptr) {
		std::cout << "[CTAnalyzer Error] FBP3D tried to project_backward_cone_fdk_from_sinograms() but FBP3D.sinograms was nullptr" << std::endl;
		ASSERT(false);
	}

	generate_blank_volume(volume_resolution_horizontal, volume_resolution_vertical, volume_resolution_horizontal);

	//ComputeProgram& kernel = *cp_ct3d_radon_transform_inverse_float_cone_equidistance;
	ComputeProgram& kernel = *cp_ct3d_radon_transform_inverse_float_cone_equidistance_sinogram_fast;

	int height = volume->get_size().y;
	for (int invocation_id = 0; invocation_id < height; invocation_id++) {
		
		kernel.update_uniform("source_sinograms", *sinograms);
		kernel.update_uniform_as_image("target_volume", *volume, 0);
		kernel.update_uniform("volume_resolution", volume->get_size());
		kernel.update_uniform("sinogram_resolution", sinograms->get_size().x, sinograms->get_size().y);
		kernel.update_uniform("sinogram_count", sinograms->get_size().z);
		kernel.update_uniform("projection_plane_width", detector_width);
		kernel.update_uniform("volume_width", volume_width);
		kernel.update_uniform("center_to_detector_distance", detector_panel_distance);
		kernel.update_uniform("center_to_source_distance", xray_source_distance);
		kernel.update_uniform("full_rotation_count", full_rotation_count);
		kernel.update_uniform("rotation_offset_radian", rotation_offset_radian);
		kernel.update_uniform("slice_y_index", invocation_id);

		kernel.dispatch(std::ceil(volume->get_size().x / 8.0f), std::ceil(volume->get_size().y / 8.0f), std::ceil(sinograms->get_size().y / 1));
	}

	clip_negatives_of_volume();
}

void FBP3D::render_sphere(glm::vec3 center, float radius, float value)
{
	center = (center / 2.0f + 0.5f) * glm::vec3(volume->get_size());
	radius = (radius / 2.0f) * (float)volume->get_size().x;

	cp_ct3d_render_sphere->update_uniform_as_image("density_volume", *volume, 0);
	cp_ct3d_render_sphere->update_uniform("volume_resolution", volume->get_size());
	cp_ct3d_render_sphere->update_uniform("sphere_center", center);
	cp_ct3d_render_sphere->update_uniform("radius", radius);
	cp_ct3d_render_sphere->update_uniform("value", glm::vec4(value, value, value, value));

	cp_ct3d_render_sphere->dispatch(std::ceil(volume->get_size().x / 4.0f), std::ceil(volume->get_size().y / 4.0f), std::ceil(volume->get_size().z / 4.0f));
}

void FBP3D::render_ellipsoid(glm::vec3 center, glm::vec3 axis_lengths, glm::vec3 rotation_angles_radian, float value)
{
	glm::vec3 volume_size = volume->get_size();
	
	center = (center / 2.0f + 0.5f) * volume_size;
	axis_lengths = (axis_lengths / 2.0f) * volume_size;

	cp_ct3d_render_ellipsoid->update_uniform_as_image("density_volume", *volume, 0);
	cp_ct3d_render_ellipsoid->update_uniform("volume_resolution", volume->get_size());
	cp_ct3d_render_ellipsoid->update_uniform("ellipsoid_center", center);
	cp_ct3d_render_ellipsoid->update_uniform("a", axis_lengths.x);
	cp_ct3d_render_ellipsoid->update_uniform("b", axis_lengths.y);
	cp_ct3d_render_ellipsoid->update_uniform("c", axis_lengths.z);
	cp_ct3d_render_ellipsoid->update_uniform("rotation_a", rotation_angles_radian.x);
	cp_ct3d_render_ellipsoid->update_uniform("rotation_b", rotation_angles_radian.y);
	cp_ct3d_render_ellipsoid->update_uniform("rotation_c", rotation_angles_radian.z);
	cp_ct3d_render_ellipsoid->update_uniform("value", glm::vec4(value, value, value, value));

	cp_ct3d_render_ellipsoid->dispatch(std::ceil(volume->get_size().x / 4.0f), std::ceil(volume->get_size().y / 4.0f), std::ceil(volume->get_size().z / 4.0f));
}

void FBP3D::blur_volume()
{
	cp_ct3d_blur_volume->update_uniform_as_image("source_density_volume", *volume, 0);
	cp_ct3d_blur_volume->update_uniform_as_image("target_density_volume", *volume, 0);
	cp_ct3d_blur_volume->update_uniform("volume_resolution", volume->get_size());

	cp_ct3d_blur_volume->dispatch(std::ceil(volume->get_size().x / 4.0f), std::ceil(volume->get_size().y / 4.0f), std::ceil(volume->get_size().z / 4.0f));
}

void FBP3D::load_volume_slice_x(int x_index, Texture2D& target_texture_out)
{
	cp_ct3d_load_x_axis->update_uniform_as_image("source_volume", *volume, 0);
	cp_ct3d_load_x_axis->update_uniform_as_image("target_texture", target_texture_out, 0);
	cp_ct3d_load_x_axis->update_uniform("volume_resolution", volume->get_size());
	cp_ct3d_load_x_axis->update_uniform("target_texture_resolution", target_texture_out.get_size());
	cp_ct3d_load_x_axis->update_uniform("slice_index", x_index);

	cp_ct3d_load_x_axis->dispatch(std::ceil(target_texture_out.get_size().x / 8.0f), std::ceil(target_texture_out.get_size().y / 8.0f), 1);
}

void FBP3D::load_volume_slice_y(int y_index, Texture2D& target_texture_out)
{
	cp_ct3d_load_y_axis->update_uniform_as_image("source_volume", *volume, 0);
	cp_ct3d_load_y_axis->update_uniform_as_image("target_texture", target_texture_out, 0);
	cp_ct3d_load_y_axis->update_uniform("volume_resolution", volume->get_size());
	cp_ct3d_load_y_axis->update_uniform("target_texture_resolution", target_texture_out.get_size());
	cp_ct3d_load_y_axis->update_uniform("slice_index", y_index);

	cp_ct3d_load_y_axis->dispatch(std::ceil(target_texture_out.get_size().x / 8.0f), std::ceil(target_texture_out.get_size().y / 8.0f), 1);
}

void FBP3D::load_volume_slice_z(int z_index, Texture2D& target_texture_out)
{
	cp_ct3d_load_z_axis->update_uniform_as_image("source_volume", *volume, 0);
	cp_ct3d_load_z_axis->update_uniform_as_image("target_texture", target_texture_out, 0);
	cp_ct3d_load_z_axis->update_uniform("volume_resolution", volume->get_size());
	cp_ct3d_load_z_axis->update_uniform("target_texture_resolution", target_texture_out.get_size());
	cp_ct3d_load_z_axis->update_uniform("slice_index", z_index);

	cp_ct3d_load_z_axis->dispatch(std::ceil(target_texture_out.get_size().x / 8.0f), std::ceil(target_texture_out.get_size().y / 8.0f), 1);
}

void FBP3D::store_volume_slice_x(Texture2D& source_texture, int x_index)
{
	cp_ct3d_store_x_axis->update_uniform_as_image("target_volume", *volume, 0);
	cp_ct3d_store_x_axis->update_uniform_as_image("source_texture", source_texture, 0);
	cp_ct3d_store_x_axis->update_uniform("volume_resolution", volume->get_size());
	cp_ct3d_store_x_axis->update_uniform("target_texture_resolution", source_texture.get_size());
	cp_ct3d_store_x_axis->update_uniform("slice_index", x_index);

	cp_ct3d_store_x_axis->dispatch(std::ceil(source_texture.get_size().x / 8.0f), std::ceil(source_texture.get_size().y / 8.0f), 1);
}

void FBP3D::store_volume_slice_y(Texture2D& source_texture, int y_index)
{
	cp_ct3d_store_y_axis->update_uniform_as_image("target_volume", *volume, 0);
	cp_ct3d_store_y_axis->update_uniform_as_image("source_texture", source_texture, 0);
	cp_ct3d_store_y_axis->update_uniform("volume_resolution", volume->get_size());
	cp_ct3d_store_y_axis->update_uniform("target_texture_resolution", source_texture.get_size());
	cp_ct3d_store_y_axis->update_uniform("slice_index", y_index);

	cp_ct3d_store_y_axis->dispatch(std::ceil(source_texture.get_size().x / 8.0f), std::ceil(source_texture.get_size().y / 8.0f), 1);

}

void FBP3D::store_volume_slice_z(Texture2D& source_texture, int z_index)
{
	cp_ct3d_store_z_axis->update_uniform_as_image("target_volume", *volume, 0);
	cp_ct3d_store_z_axis->update_uniform_as_image("source_texture", source_texture, 0);
	cp_ct3d_store_z_axis->update_uniform("volume_resolution", volume->get_size());
	cp_ct3d_store_z_axis->update_uniform("target_texture_resolution", source_texture.get_size());
	cp_ct3d_store_z_axis->update_uniform("slice_index", z_index);

	cp_ct3d_store_z_axis->dispatch(std::ceil(source_texture.get_size().x / 8.0f), std::ceil(source_texture.get_size().y / 8.0f), 1);
}

void FBP3D::load_projection(int projection_index, Texture2D& target_texture_out)
{
	cp_ct3d_projection_load->update_uniform_as_image("source_projections", *projections, 0);
	cp_ct3d_projection_load->update_uniform_as_image("target_texture", target_texture_out, 0);
	cp_ct3d_projection_load->update_uniform("projection_resolution", projections->get_size());
	cp_ct3d_projection_load->update_uniform("target_texture_resolution", target_texture_out.get_size());
	cp_ct3d_projection_load->update_uniform("projection_index", projection_index);

	cp_ct3d_projection_load->dispatch(std::ceil(target_texture_out.get_size().x / 8.0f), std::ceil(target_texture_out.get_size().y / 8.0f), 1);

}

void FBP3D::store_projection(Texture2D& source_texture, int projection_index)
{
	cp_ct3d_projection_store->update_uniform_as_image("source_texture", source_texture, 0);
	cp_ct3d_projection_store->update_uniform_as_image("target_projections", *projections, 0);
	cp_ct3d_projection_store->update_uniform("projection_resolution", projections->get_size());
	cp_ct3d_projection_store->update_uniform("source_texture_resolution", source_texture.get_size());
	cp_ct3d_projection_store->update_uniform("projection_index", projection_index);

	cp_ct3d_projection_store->dispatch(std::ceil(projections->get_size().x / 8.0f), std::ceil(projections->get_size().y / 8.0f), 1);

}

void FBP3D::load_sinogram(int y_coord, Texture2D& target_texture_out)
{
	ComputeProgram& kernel = *cp_ct3d_sinograms_load_z_axis;
	kernel.update_uniform_as_image("source_texture", *sinograms, 0);
	kernel.update_uniform_as_image("target_texture", target_texture_out, 0);
	kernel.update_uniform("array_texture_resolution", sinograms->get_size());
	kernel.update_uniform("target_texture_resolution", target_texture_out.get_size());
	kernel.update_uniform("slice_index", y_coord);

	kernel.dispatch(std::ceil(target_texture_out.get_size().x / 8.0f), std::ceil(target_texture_out.get_size().y / 8.0f), 1);
}

void FBP3D::store_sinogram(Texture2D& soruce_texture, int y_coord)
{
	ComputeProgram& kernel = *cp_ct3d_sinograms_store_z_axis;
	kernel.update_uniform_as_image("target_texture", *sinograms, 0);
	kernel.update_uniform_as_image("source_texture", soruce_texture, 0);
	kernel.update_uniform("array_texture_resolution", sinograms->get_size());
	kernel.update_uniform("target_texture_resolution", soruce_texture.get_size());
	kernel.update_uniform("slice_index", y_coord);

	kernel.dispatch(std::ceil(sinograms->get_size().x / 8.0f), std::ceil(sinograms->get_size().y / 8.0f), 1);

}

void FBP3D::load_projection_from_sinograms(int projection_index, Texture2D& target_texture_out)
{
	ComputeProgram& kernel = *cp_ct3d_sinograms_load_z_axis;
	kernel.update_uniform_as_image("source_texture", *sinograms, 0);
	kernel.update_uniform_as_image("target_texture", target_texture_out, 0);
	kernel.update_uniform("array_texture_resolution", sinograms->get_size());
	kernel.update_uniform("target_texture_resolution", target_texture_out.get_size());
	kernel.update_uniform("projection_index", projection_index);

	kernel.dispatch(std::ceil(target_texture_out.get_size().x / 8.0f), std::ceil(target_texture_out.get_size().y / 8.0f), 1);

}

void FBP3D::store_projection_to_sinograms(Texture2D& soruce_texture, int projection_index)
{
	ComputeProgram& kernel = *cp_ct3d_sinograms_store_y_axis;
	kernel.update_uniform_as_image("target_texture", *sinograms, 0);
	kernel.update_uniform_as_image("source_texture", soruce_texture, 0);
	kernel.update_uniform("array_texture_resolution", sinograms->get_size());
	kernel.update_uniform("target_texture_resolution", soruce_texture.get_size());
	kernel.update_uniform("projection_index", projection_index);

	kernel.dispatch(std::ceil(sinograms->get_size().x / 8.0f), std::ceil(sinograms->get_size().y / 8.0f), 1);
}

void FBP3D::compute_max_value_of_volume(Texture1D& target_int_texture)
{
	ComputeProgram& kernel = *cp_ct3d_compute_max_to_int_texture;

	kernel.update_uniform_as_image("source_volume", *volume, 0);
	kernel.update_uniform_as_image("target_max_texture", target_int_texture, 0);
	kernel.update_uniform("source_volume_resolution", volume->get_size());

	kernel.dispatch(std::ceil(volume->get_size().x / 8.0), std::ceil(volume->get_size().y / 8.0), std::ceil(volume->get_size().z / 1.0));
}

void FBP3D::compute_min_value_of_volume(Texture1D& target_int_texture)
{
	ComputeProgram& kernel = *cp_ct3d_compute_min_to_int_texture;

	kernel.update_uniform_as_image("source_volume", *volume, 0);
	kernel.update_uniform_as_image("target_min_texture", target_int_texture, 0);
	kernel.update_uniform("source_volume_resolution", volume->get_size());

	kernel.dispatch(std::ceil(volume->get_size().x / 8.0), std::ceil(volume->get_size().y / 8.0), std::ceil(volume->get_size().z / 1.0));
}

void FBP3D::normalize_histogram()
{
	std::shared_ptr<Texture1D> min_int_texture = std::make_shared<Texture1D>(1, Texture1D::ColorTextureFormat::R32UI, 1, 0.0f);
	min_int_texture->is_bindless = false;

	std::shared_ptr<Texture1D> max_int_texture = std::make_shared<Texture1D>(1, Texture1D::ColorTextureFormat::R32UI, 1, 0.0f);
	max_int_texture->is_bindless = false;

	compute_max_value_of_volume(*max_int_texture);
	compute_min_value_of_volume(*min_int_texture);

	normalize_histogram(*min_int_texture, *max_int_texture);
}

void FBP3D::normalize_histogram(Texture1D& source_min_int_texture, Texture1D& source_max_int_texture)
{
	ComputeProgram& kernel = *cp_ct3d_histogram_normalize;

	kernel.update_uniform_as_image("source_volume", *volume, 0);
	kernel.update_uniform_as_image("source_min_texture", source_min_int_texture, 0);
	kernel.update_uniform_as_image("source_max_texture", source_max_int_texture, 0);
	kernel.update_uniform("source_volume_resolution", volume->get_size());

	kernel.dispatch(std::ceil(volume->get_size().x / 8.0), std::ceil(volume->get_size().y / 8.0), std::ceil(volume->get_size().z / 1.0));
}
