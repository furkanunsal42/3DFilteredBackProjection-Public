#include "GraphicsCortex.h"
#include "Data/CTData.h"

int main() {

	glm::ivec3 volume_dimentions(512, 512, 512);
	int projection_count = 512;
	int window_width = 1024;

	Frame frame(window_width, window_width, "CTAnalyzer", 0, 1, true, true, false, Frame::CallbackLevel::NOTIFICATION, false);
	Scene scene(frame);
	scene.camera->fov = 90;
	scene.camera->max_distance = 1000;

	std::shared_ptr<CTData> data = std::make_shared<CTData>();
	data->generate_shepplogan(volume_dimentions.x, volume_dimentions.y, volume_dimentions.z);

	data->project_forward_parallel(1, projection_count, 0);
	data->apply_filter_to_projections(FBP::FilterType::COSINE);
	data->project_backward_parallel(1, volume_dimentions.x, 0);
	data->project_forward_parallel(1, projection_count, 0);

	std::shared_ptr<Texture2D> slice = std::make_shared<Texture2D>(volume_dimentions.x, volume_dimentions.y, Texture2D::ColorTextureFormat::R32F, 1, 0, 0);
	slice->is_bindless = false;
	std::shared_ptr<Texture2D> slice_complex = std::make_shared<Texture2D>(volume_dimentions.x, volume_dimentions.y, Texture2D::ColorTextureFormat::RG32F, 1, 0, 0);
	slice_complex->is_bindless = false;
	std::shared_ptr<Texture2D> slice_white = std::make_shared<Texture2D>(volume_dimentions.x, volume_dimentions.y, Texture2D::ColorTextureFormat::RGBA32F, 1, 0, 0);
	slice_white->is_bindless = false;

	std::shared_ptr<Framebuffer> framebuffer = std::make_shared<Framebuffer>();

	int i = 0;
	while (frame.is_running()) {
		double delta_time = frame.handle_window();
		frame.clear_window();
		frame.display_performance(180);

		i++;
		i = i % volume_dimentions.z;

		data->load_volume_slice_y(i, *slice);
		//data->load_projection(i, *slice);

		data->fbp_solver->_texture_blit_float1_to_float4(*slice, *slice_white);

		framebuffer->attach_color(0, slice_white);
		framebuffer->set_read_buffer(0);
		framebuffer->blit_to_screen(0, 0, volume_dimentions.x, volume_dimentions.y, 0, 0, window_width, window_width, Framebuffer::Channel::COLOR, Framebuffer::Filter::LINEAR);
	}
}