#pragma once
#include <string>
#include <unordered_map>
#include <memory>

#include "Texture1D.h"
#include "Texture2D.h"
#include "Texture2DArray.h"
#include "Texture3D.h"

#include "ComputeProgram.h"

#include "FastFourierTransformGPU/FastFourierTransformGPU.h"
#include "FilteredBackProjection/FBP2D.h"

class FBP3D {
public:
	
	FBP3D(const FBP3D& other) = delete;
	FBP3D();

	void generate_blank_volume(int x_size, int y_size, int z_size);
	void generate_blank_projections(int x_size, int y_size, int projection_amount);
	void generate_blank_sinograms(int x_size, int y_size, int projection_amount);

	void render_sphere(glm::vec3 center, float radius, float value);
	void render_ellipsoid(glm::vec3 center, glm::vec3 axis_lengths, glm::vec3 rotation_angles_radian, float value);
	void blur_volume();
	void generate_shepplogan(int x_size, int y_size, int z_size);
	
	void read_volume(const std::string& upper_diretory, int slice_width, int slice_height, int byte_per_channel, int slice_count);
	void read_projections(const std::string& upper_diretory, int file_width, int file_height, int file_channel_count, int byte_per_channel, int target_projection_width, int target_projection_height, int projection_count);
	void read_sinograms();
	void read_projections_as_sinograms(const std::string& upper_diretory, int file_width, int file_height, int file_channel_count, int byte_per_channel, int target_projection_width, int target_projection_height, int projection_count);

	void write_volume(std::string upper_directory);
	void write_projections(const std::string& upper_directory);
	void write_sinograms();

	void free_volume();
	void free_projections();
	void free_sinograms();

	void project_forward_parallel_to_projections(float full_rotation_count, int projection_count, float rotation_offset_radian = 0.0f);
	void project_forward_cone_to_projections(float xray_source_distance, float detector_panel_distance, float detector_width, float volume_width, float full_rotation_count, int projection_count, float rotation_offset_radian = 0.0f);
	void project_forward_parallel_to_sinogram(float full_rotation_count, int projection_count, float rotation_offset_radian = 0.0f);
	void project_forward_cone_to_sinogram(float xray_source_distance, float detector_panel_distance, float detector_width, float volume_width, float full_rotation_count, int projection_count, float rotation_offset_radian = 0.0f);

	void log_normalize_projections(float base_value);
	void apply_fdk_weights_to_projections(float xray_source_distance, float detector_panel_distance, float detector_width);
	void apply_filter_to_projections(FBP2D::FilterType filter_type);
	void log_normalize_sinograms(float base_value);
	void apply_fdk_weights_to_sinograms(float xray_source_distance, float detector_panel_distance, float detector_width);
	void apply_filter_to_sinograms(FBP2D::FilterType filter_type);
	void clip_negatives_of_volume();
	
	void project_backward_parallel_from_projections(float full_rotation_count, int volume_resolution_horizontal, int volume_resolution_vertical, float rotation_offset_radian = 0.0f);
	void project_backward_cone_fdk_from_projections(float xray_source_distance, float detector_panel_distance, float detector_width, float volume_width, float full_rotation_count, int volume_resolution_horizontal, int volume_resolution_vertical, float rotation_offset_radian = 0.0f);
	void project_backward_parallel_from_sinograms(float full_rotation_count, int volume_resolution_horizontal, int volume_resolution_vertical, float rotation_offset_radian = 0.0f);
	void project_backward_cone_fdk_from_sinograms(float xray_source_distance, float detector_panel_distance, float detector_width, float volume_width, float full_rotation_count, int volume_resolution_horizontal, int volume_resolution_vertical, float rotation_offset_radian = 0.0f);

	void load_volume_slice_x(int x_index, Texture2D& target_texture_out);
	void load_volume_slice_y(int y_index, Texture2D& target_texture_out);
	void load_volume_slice_z(int z_index, Texture2D& target_texture_out);

	void store_volume_slice_x(Texture2D& source_texture, int x_index);
	void store_volume_slice_y(Texture2D& source_texture, int y_index);
	void store_volume_slice_z(Texture2D& source_texture, int z_index);

	void load_projection(int projection_index, Texture2D& target_texture_out);
	void store_projection(Texture2D& source_texture, int projection_index);

	void load_sinogram(int projection_index, Texture2D& target_texture_out);
	void store_sinogram(Texture2D& source_texture, int projection_index);

	void load_projection_from_sinograms(int projection_index, Texture2D& target_texture_out);
	void store_projection_to_sinograms(Texture2D& soruce_texture, int projection_index);

	void compute_max_value_of_volume(Texture1D& target_int_texture);
	void compute_min_value_of_volume(Texture1D& target_int_texture);
	void normalize_histogram();
	void normalize_histogram(Texture1D& source_min_int_texture, Texture1D& source_max_int_texture);

	std::shared_ptr<FBP2D> fbp_solver = std::make_shared<FBP2D>();
	std::shared_ptr<FFFT> fft_solver = std::make_shared<FFFT>();

	std::shared_ptr<Texture3D> volume = nullptr;
	std::shared_ptr<Texture2DArray> projections = nullptr;
	std::shared_ptr<Texture3D> sinograms = nullptr;

	std::shared_ptr<ComputeProgram> cp_ct3d_render_sphere		= std::make_shared<ComputeProgram>
		(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/ct3d_render_sphere.comp"));
	std::shared_ptr<ComputeProgram> cp_ct3d_render_ellipsoid	= std::make_shared<ComputeProgram>
		(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/ct3d_render_ellipsoid.comp"));
	std::shared_ptr<ComputeProgram> cp_ct3d_blur_volume			= std::make_shared<ComputeProgram>
		(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/ct3d_blur_volume.comp"));
	std::shared_ptr<ComputeProgram> cp_ct3d_log_normalize		= std::make_shared<ComputeProgram>
		(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/ct3d_log_normalize.comp"));
	std::shared_ptr<ComputeProgram> cp_ct3d_fdk_weight			= std::make_shared<ComputeProgram>
		(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/ct3d_fdk_weight.comp"));
	std::shared_ptr<ComputeProgram> cp_ct3d_clip_negative		= std::make_shared<ComputeProgram>
		(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/ct3d_clip_negative.comp"));

	std::shared_ptr<ComputeProgram> cp_ct3d_compute_min_to_int_texture	= std::make_shared<ComputeProgram>
		(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/ct3d_compute_min_to_int_texture.comp"));
	std::shared_ptr<ComputeProgram> cp_ct3d_compute_max_to_int_texture	= std::make_shared<ComputeProgram>
		(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/ct3d_compute_max_to_int_texture.comp"));
	std::shared_ptr<ComputeProgram> cp_ct3d_histogram_normalize			= std::make_shared<ComputeProgram>
		(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/ct3d_histogram_normalize.comp"));

	std::shared_ptr<ComputeProgram> cp_ct3d_load_x_axis			= std::make_shared<ComputeProgram>
		(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/ct3d_load_x_axis.comp"));
	std::shared_ptr<ComputeProgram> cp_ct3d_load_y_axis			= std::make_shared<ComputeProgram>
		(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/ct3d_load_y_axis.comp"));
	std::shared_ptr<ComputeProgram> cp_ct3d_load_z_axis			= std::make_shared<ComputeProgram>
		(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/ct3d_load_z_axis.comp"));
	std::shared_ptr<ComputeProgram> cp_ct3d_store_x_axis		= std::make_shared<ComputeProgram>
		(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/ct3d_store_x_axis.comp"));
	std::shared_ptr<ComputeProgram> cp_ct3d_store_y_axis		= std::make_shared<ComputeProgram>
		(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/ct3d_store_y_axis.comp"));
	std::shared_ptr<ComputeProgram> cp_ct3d_store_z_axis		= std::make_shared<ComputeProgram>
		(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/ct3d_store_z_axis.comp"));
	std::shared_ptr<ComputeProgram> cp_ct3d_projection_load		= std::make_shared<ComputeProgram>
		(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/ct3d_projection_load.comp"));
	std::shared_ptr<ComputeProgram> cp_ct3d_projection_store	= std::make_shared<ComputeProgram>
		(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/ct3d_projection_store.comp"));
	std::shared_ptr<ComputeProgram> cp_ct3d_sinograms_store_z_axis = std::make_shared<ComputeProgram>
		(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/ct3d_sinograms_store_z_axis.comp"));
	std::shared_ptr<ComputeProgram> cp_ct3d_sinograms_load_z_axis = std::make_shared<ComputeProgram>
		(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/ct3d_sinograms_load_z_axis.comp"));
	std::shared_ptr<ComputeProgram> cp_ct3d_sinograms_store_y_axis = std::make_shared<ComputeProgram>
		(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/ct3d_sinograms_store_y_axis.comp"));
	std::shared_ptr<ComputeProgram> cp_ct3d_sinograms_load_y_axis = std::make_shared<ComputeProgram>
		(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/ct3d_sinograms_load_y_axis.comp"));

	std::shared_ptr<ComputeProgram> cp_ct3d_radon_transform_float_parallel							= std::make_shared<ComputeProgram>
		(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/ct3d_radon_transform_float_parallel.comp"));
	std::shared_ptr<ComputeProgram> cp_ct3d_radon_transform_inverse_float_parallel					= std::make_shared<ComputeProgram>
		(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/ct3d_radon_transform_inverse_float_parallel.comp"));
	std::shared_ptr<ComputeProgram> cp_ct3d_radon_transform_float_cone_equidistance					= std::make_shared<ComputeProgram>
		(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/ct3d_radon_transform_float_cone_equidistance.comp"));
	std::shared_ptr<ComputeProgram> cp_ct3d_radon_transform_inverse_float_cone_equidistance			= std::make_shared<ComputeProgram>
		(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/ct3d_radon_transform_inverse_float_cone_equidistance.comp"));
	std::shared_ptr<ComputeProgram> cp_ct3d_radon_transform_inverse_float_cone_equidistance_fast	= std::make_shared<ComputeProgram>
		(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/ct3d_radon_transform_inverse_float_cone_equidistance_fast.comp"));

	std::shared_ptr<ComputeProgram> cp_ct3d_radon_transform_float_parallel_sinogram							= std::make_shared<ComputeProgram>
		(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/ct3d_radon_transform_float_parallel_sinogram.comp"));
	std::shared_ptr<ComputeProgram> cp_ct3d_radon_transform_float_cone_equidistance_sinogram				= std::make_shared<ComputeProgram>
		(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/ct3d_radon_transform_float_cone_equidistance_sinogram.comp"));
	std::shared_ptr<ComputeProgram> cp_ct3d_radon_transform_inverse_float_parallel_sinogram					= std::make_shared<ComputeProgram>
		(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/ct3d_radon_transform_inverse_float_parallel_sinogram.comp"));
	std::shared_ptr<ComputeProgram> cp_ct3d_radon_transform_inverse_float_cone_equidistance_sinogram_fast	= std::make_shared<ComputeProgram>
		(Shader("../CTAnalyzer/Source/GLSL/Compute/FBP/ct3d_radon_transform_inverse_float_cone_equidistance_sinogram_fast.comp"));

	const Texture2D::ColorTextureFormat _volume_internal_format				= Texture2D::ColorTextureFormat::R32F;
	const Texture2D::ColorTextureFormat _projection_internal_format			= Texture2D::ColorTextureFormat::R32F;
	const Texture2D::ColorTextureFormat _projection_complex_internal_format = Texture2D::ColorTextureFormat::RG32F;

private:
	const int _slice_index_number_length = 6;
	
	const int _volume_slices_batch_height = 512;
	const int _slice_array_depth_per_batch = 32;

	int _slice_width;
	int _slice_height;
	int _slice_bytes_per_channel;
	int _projection_count;
};