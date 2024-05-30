#include "FastFourierTransformGPU/FastFourierTransformGPU.h"

void FFFT::fft_radix2(Texture2D& complex_texture)
{
	reverse_bit_swap(complex_texture);

	int size = complex_texture.get_size().x;
	int log2_size = floor_log2(size);

	if (1U << log2_size != size) {
		std::cout << "[CTAnalyzer Error] FFFT::fft_radix2() is called but texture's width is not a power of 2" << std::endl;
		ASSERT(false);
	}

	std::shared_ptr<Texture2D> complex_texture_copy = complex_texture.create_texture_with_same_parameters();
	blit_texture_complex_to_complex(complex_texture, *complex_texture_copy);

	for (int i = 0; i < log2_size; i++) {
		// swap the read and write buffers every iteration to avoid rece-condition and unnecessary copying
		Texture2D& texture_to_read  = i % 2 == 0 ? complex_texture : *complex_texture_copy;
		Texture2D& texture_to_write = i % 2 == 0 ? *complex_texture_copy : complex_texture;

		fft_single_step(texture_to_read, texture_to_write, i);
	}

	if ((log2_size - 1) % 2 == 0)
		blit_texture_complex_to_complex(*complex_texture_copy, complex_texture);
}

void FFFT::inverse_fft_radix2(Texture2D& complex_texture)
{
	conjugate_complex_texture(complex_texture);
	fft_radix2(complex_texture);
	conjugate_complex_texture(complex_texture);
	divide_complex_texture(complex_texture, (float)complex_texture.get_size().x);
}

void FFFT::fft(Texture2D& source_complex_texture, Texture2D& target_fourier_texture)
{
	int size = source_complex_texture.get_size().x;
	int logsize = 1 << ceil_log2(size);

	if (target_fourier_texture.get_size().x != logsize) {
		std::cout << "[CTAnalyzer Error] FFFT::fft(Texture2D& source_complex_texture, Texture2D& target_fourier_texture) was called but target_fourier_texture's x dimention is now radix-2, it had to be " << logsize << std::endl;
		ASSERT(false);
	}

	target_fourier_texture.clear(glm::vec2(0, 0), 0);
	blit_texture_complex_to_complex(source_complex_texture, target_fourier_texture);
	
	fft_radix2(target_fourier_texture);
}

std::shared_ptr<Texture2D> FFFT::fft(Texture2D& source_complex_texture)
{
	int size = source_complex_texture.get_size().x;
	int logsize = 1 << ceil_log2(size);

	std::shared_ptr<Texture2D> target_fourier_texture = created_padded_complex_texture(source_complex_texture);

	fft(source_complex_texture, *target_fourier_texture);

	return target_fourier_texture;
}

void FFFT::inverse_fft(Texture2D& source_fourier_texture, Texture2D& target_texture)
{
	if (source_fourier_texture.get_size().x == target_texture.get_size().x) {
		blit_texture_complex_to_complex(source_fourier_texture, target_texture);

		inverse_fft_radix2(target_texture);
	}
	else {
		std::shared_ptr<Texture2D> source_fourier_copy = source_fourier_texture.create_texture_with_same_parameters();
		blit_texture_complex_to_complex(source_fourier_texture, *source_fourier_copy);

		inverse_fft_radix2(*source_fourier_copy);

		blit_texture_complex_to_complex(*source_fourier_copy, target_texture);
	}
}

std::shared_ptr<Texture2D> FFFT::inverse_fft(Texture2D& source_fourier_texture, int texture_size_before_padding)
{
	std::shared_ptr<Texture2D> target_fourier_texture = std::make_shared<Texture2D>(
		texture_size_before_padding, source_fourier_texture.get_size().y, Texture2D::ColorTextureFormat::RG32F, 1, 0, 0);

	inverse_fft(source_fourier_texture, *target_fourier_texture);

	return target_fourier_texture;
}

void FFFT::fft_shift(Texture2D& source_complex_texture, Texture2D& target_complex_texture)
{
	fft_shift(source_complex_texture, target_complex_texture, source_complex_texture.get_size().x / 2);
}

void FFFT::fft_shift(Texture2D& complex_texture)
{
	fft_shift(complex_texture, complex_texture.get_size().x / 2);
}

void FFFT::fft_shift(Texture2D& source_complex_texture, Texture2D& target_complex_texture, int shift_amount)
{
	cp_shift->update_uniform_as_image("source_complex_texture", source_complex_texture, 0);
	cp_shift->update_uniform_as_image("target_complex_texture", target_complex_texture, 0);
	cp_shift->update_uniform("texture_resolution", source_complex_texture.get_size());
	cp_shift->update_uniform("shift_amount", shift_amount);

	cp_shift->dispatch(std::ceil(source_complex_texture.get_size().x / 8.0f), std::ceil(source_complex_texture.get_size().y / 8.0f), 1);
}

void FFFT::fft_shift(Texture2D& complex_texture, int shift_amount)
{
	std::shared_ptr<Texture2D> copy_texture = complex_texture.create_texture_with_same_parameters();
	
	fft_shift(complex_texture, *copy_texture, shift_amount);
	blit_texture_complex_to_complex(*copy_texture, complex_texture);
}

void FFFT::inverse_fft_shift(Texture2D& source_complex_texture, Texture2D& target_complex_texture)
{
	inverse_fft_shift(source_complex_texture, target_complex_texture, source_complex_texture.get_size().x / 2);
}

void FFFT::inverse_fft_shift(Texture2D& complex_texture)
{
	inverse_fft_shift(complex_texture, complex_texture.get_size().x / 2);
}

void FFFT::inverse_fft_shift(Texture2D& source_complex_texture, Texture2D& target_complex_texture, int shift_amount)
{
	fft_shift(source_complex_texture, target_complex_texture, -shift_amount);
}

void FFFT::inverse_fft_shift(Texture2D& complex_texture, int shift_amount)
{
	std::shared_ptr<Texture2D> copy_texture = complex_texture.create_texture_with_same_parameters();

	inverse_fft_shift(complex_texture, *copy_texture, shift_amount);
	blit_texture_complex_to_complex(*copy_texture, complex_texture);
}

void FFFT::blit_texture_complex_to_complex(Texture2D& source_complex_texture, Texture2D& target_complex_texture)
{
	glm::ivec2 size = glm::min(source_complex_texture.get_size(), target_complex_texture.get_size());

	cp_blit_texture_complex_to_complex->update_uniform_as_image("source_texture", source_complex_texture, 0);
	cp_blit_texture_complex_to_complex->update_uniform_as_image("target_texture", target_complex_texture, 0);
	cp_blit_texture_complex_to_complex->update_uniform("texture_resolution", size);

	cp_blit_texture_complex_to_complex->dispatch(std::ceil(target_complex_texture.get_size().x / 8.0f), std::ceil(target_complex_texture.get_size().y / 8.0f), 1);
}

void FFFT::blit_texture_real_to_complex(Texture2D& source_real_texture, Texture2D& target_complex_texture)
{
	glm::ivec2 size = glm::min(source_real_texture.get_size(), target_complex_texture.get_size());
	
	cp_blit_texture_real_to_complex->update_uniform_as_image("source_texture", source_real_texture, 0);
	cp_blit_texture_real_to_complex->update_uniform_as_image("target_texture", target_complex_texture, 0);
	cp_blit_texture_real_to_complex->update_uniform("texture_resolution", size);

	cp_blit_texture_real_to_complex->dispatch(std::ceil(target_complex_texture.get_size().x / 8.0f), std::ceil(target_complex_texture.get_size().y / 8.0f), 1);
}

void FFFT::blit_texture_complex_to_real_r(Texture2D& source_complex_texture, Texture2D& target_real_texture)
{
	glm::ivec2 size = glm::min(source_complex_texture.get_size(), target_real_texture.get_size());
	
	cp_blit_texture_complex_to_real_r->update_uniform_as_image("source_texture", source_complex_texture, 0);
	cp_blit_texture_complex_to_real_r->update_uniform_as_image("target_texture", target_real_texture, 0);
	cp_blit_texture_complex_to_real_r->update_uniform("texture_resolution", size);

	cp_blit_texture_complex_to_real_r->dispatch(std::ceil(target_real_texture.get_size().x / 8.0f), std::ceil(target_real_texture.get_size().y / 8.0f), 1);
}

void FFFT::blit_texture_complex_to_real_i(Texture2D& source_complex_texture, Texture2D& target_real_texture)
{
	glm::ivec2 size = glm::min(source_complex_texture.get_size(), target_real_texture.get_size());
	
	cp_blit_texture_complex_to_real_i->update_uniform_as_image("source_texture", source_complex_texture, 0);
	cp_blit_texture_complex_to_real_i->update_uniform_as_image("target_texture", target_real_texture, 0);
	cp_blit_texture_complex_to_real_i->update_uniform("texture_resolution", size);

	cp_blit_texture_complex_to_real_i->dispatch(std::ceil(target_real_texture.get_size().x / 8.0f), std::ceil(target_real_texture.get_size().y / 8.0f), 1);
}

std::shared_ptr<Texture2D> FFFT::create_complex_texture_from_real(Texture2D& real_texture)
{
	std::shared_ptr<Texture2D> complex_texture = std::make_shared<Texture2D>(real_texture.get_size().x, real_texture.get_size().y, Texture2D::ColorTextureFormat::RG32F, 1, 0, 0);
	blit_texture_real_to_complex(real_texture, *complex_texture);
	
	return complex_texture;
}

std::shared_ptr<Texture2D> FFFT::create_real_texture_from_complex_r(Texture2D& complex_texture)
{
	std::shared_ptr<Texture2D> real_texture = std::make_shared<Texture2D>(complex_texture.get_size().x, complex_texture.get_size().y, Texture2D::ColorTextureFormat::R32F, 1, 0, 0);
	blit_texture_complex_to_real_r(complex_texture, *real_texture);

	return real_texture;
}

std::shared_ptr<Texture2D> FFFT::create_real_texture_from_complex_i(Texture2D& complex_texture)
{
	std::shared_ptr<Texture2D> real_texture = std::make_shared<Texture2D>(complex_texture.get_size().x, complex_texture.get_size().y, Texture2D::ColorTextureFormat::R32F, 1, 0, 0);
	blit_texture_complex_to_real_i(complex_texture, *real_texture);

	return real_texture;
}

void FFFT::conjugate_complex_texture(Texture2D& texture)
{
	cp_conjugate_texture->update_uniform_as_image("complex_texture", texture, 0);
	cp_conjugate_texture->update_uniform("texture_resolution", texture.get_size());

	cp_conjugate_texture->dispatch(std::ceil(texture.get_size().x / 8.0f), std::ceil(texture.get_size().y / 8.0f), 1);
}

void FFFT::multiply_complex_texture(Texture2D& texture, float multiplier)
{
	cp_multiply_complex_texture->update_uniform_as_image("complex_texture", texture, 0);
	cp_multiply_complex_texture->update_uniform("multiplier", multiplier);
	cp_multiply_complex_texture->update_uniform("texture_resolution", texture.get_size());

	cp_multiply_complex_texture->dispatch(std::ceil(texture.get_size().x / 8.0f), std::ceil(texture.get_size().y / 8.0f), 1);
}

void FFFT::divide_complex_texture(Texture2D& texture, float divisor)
{
	cp_divide_complex_texture->update_uniform_as_image("complex_texture", texture, 0);
	cp_divide_complex_texture->update_uniform("divisor", divisor);
	cp_divide_complex_texture->update_uniform("texture_resolution", texture.get_size());

	cp_divide_complex_texture->dispatch(std::ceil(texture.get_size().x / 8.0f), std::ceil(texture.get_size().y / 8.0f), 1);
}

void FFFT::reverse_bit_swap(Texture2D& texture)
{
	cp_reverse_bit_swap->update_uniform_as_image("target_image", texture, 0);
	cp_reverse_bit_swap->update_uniform("image_resolution", texture.get_size());
	cp_reverse_bit_swap->update_uniform("ceil_log2_size", floor_log2(texture.get_size().x));
	cp_reverse_bit_swap->dispatch(std::ceil(texture.get_size().x / 8.0f), std::ceil(texture.get_size().y / 8.0f), 1);
}

void FFFT::fft_single_step(Texture2D& read_texture, Texture2D& write_texture, int step_index)
{
	int size = 2 << step_index;

	cp_fft_single_pass_complex_complex->update_uniform_as_image("fft_read_texture", read_texture, 0);
	cp_fft_single_pass_complex_complex->update_uniform_as_image("fft_write_texture", write_texture, 0);
	cp_fft_single_pass_complex_complex->update_uniform("fft_texture_resolution", read_texture.get_size());
	cp_fft_single_pass_complex_complex->update_uniform("size", size);

	cp_fft_single_pass_complex_complex->dispatch(std::ceil(read_texture.get_size().x / 8.0f), std::ceil(read_texture.get_size().y / 8.0f), 1);
}

std::shared_ptr<Texture2D> FFFT::created_padded_complex_texture(Texture2D& source_complex_texture)
{
	int padded_width = 1 << ceil_log2(source_complex_texture.get_size().x);
	std::shared_ptr<Texture2D> padded_texture = std::make_shared<Texture2D>(padded_width, source_complex_texture.get_size().y, Texture2D::ColorTextureFormat::RG32F, 1, 0, 0);
	padded_texture->is_bindless = false;

	blit_texture_complex_to_complex(source_complex_texture, *padded_texture);

	return padded_texture;
}

int FFFT::floor_log2(int n)
{
	int log2_size = 0;
	for (size_t temp = n; temp > 1U; temp >>= 1)
		log2_size++;
	return log2_size;
}

int FFFT::ceil_log2(int n)
{
	return floor_log2(2 * n - 1);
}
