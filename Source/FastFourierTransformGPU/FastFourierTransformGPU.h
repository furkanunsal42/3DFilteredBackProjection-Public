#pragma once

#include <memory>

#include "ComputeProgram.h"
#include "Texture2D.h"

// Fucking Fast Fourier Transform
// parallel fast fourier transform on GPU
class FFFT {
public:

	void fft_radix2(Texture2D& complex_texture);
	void inverse_fft_radix2(Texture2D& complex_texture);

	void fft(Texture2D& source_complex_texture, Texture2D& target_fourier_texture);
	std::shared_ptr<Texture2D> fft(Texture2D& source_complex_texture);
	
	void inverse_fft(Texture2D& source_fourier_texture, Texture2D& target_texture);
	std::shared_ptr<Texture2D> inverse_fft(Texture2D& source_fourier_texture, int texture_size_before_padding);
	
	void fft_shift(Texture2D& source_complex_texture, Texture2D& target_complex_texture);
	void fft_shift(Texture2D& complex_texture);
	void fft_shift(Texture2D& source_complex_texture, Texture2D& target_complex_texture, int shift_amount);
	void fft_shift(Texture2D& complex_texture, int shift_amount);
	void inverse_fft_shift(Texture2D& source_complex_texture, Texture2D& target_complex_texture);
	void inverse_fft_shift(Texture2D& complex_texture);
	void inverse_fft_shift(Texture2D& source_complex_texture, Texture2D& target_complex_texture, int shift_amount);
	void inverse_fft_shift(Texture2D& complex_texture, int shift_amount);

	// helper functions
	void reverse_bit_swap(Texture2D& texture);
	void fft_single_step(Texture2D& read_texture, Texture2D& write_texture, int step_index);

	std::shared_ptr<Texture2D> created_padded_complex_texture(Texture2D& source_complex_texture);

	void blit_texture_complex_to_complex(Texture2D& source_complex_texture, Texture2D& target_complex_texture);
	void blit_texture_real_to_complex(Texture2D& source_real_texture, Texture2D& target_complex_texture);
	void blit_texture_complex_to_real_r(Texture2D& source_complex_texture, Texture2D& target_real_texture);
	void blit_texture_complex_to_real_i(Texture2D& source_complex_texture, Texture2D& target_real_texture);
	
	std::shared_ptr<Texture2D> create_complex_texture_from_real(Texture2D& real_texture);
	std::shared_ptr<Texture2D> create_real_texture_from_complex_r(Texture2D& complex_texture);
	std::shared_ptr<Texture2D> create_real_texture_from_complex_i(Texture2D& complex_texture);

	void conjugate_complex_texture(Texture2D& texture);
	void multiply_complex_texture(Texture2D& texture, float multiplier);
	void divide_complex_texture(Texture2D& texture, float divisor);

	std::shared_ptr<ComputeProgram> cp_reverse_bit_swap					= std::make_shared<ComputeProgram>(Shader("../CTAnalyzer/Source/GLSL/Compute/FFT/reverse_bit_swap.comp"));
	std::shared_ptr<ComputeProgram> cp_fft_single_pass_complex_complex	= std::make_shared<ComputeProgram>(Shader("../CTAnalyzer/Source/GLSL/Compute/FFT/fft_single_step_complex_complex.comp"));
	std::shared_ptr<ComputeProgram> cp_shift							= std::make_shared<ComputeProgram>(Shader("../CTAnalyzer/Source/GLSL/Compute/FFT/fft_shift.comp"));

	std::shared_ptr<ComputeProgram> cp_blit_texture_complex_to_complex	= std::make_shared<ComputeProgram>(Shader("../CTAnalyzer/Source/GLSL/Compute/FFT/blit_texture_complex_to_complex.comp"));
	std::shared_ptr<ComputeProgram> cp_blit_texture_real_to_complex		= std::make_shared<ComputeProgram>(Shader("../CTAnalyzer/Source/GLSL/Compute/FFT/blit_texture_real_to_complex.comp"));
	std::shared_ptr<ComputeProgram> cp_blit_texture_complex_to_real_r	= std::make_shared<ComputeProgram>(Shader("../CTAnalyzer/Source/GLSL/Compute/FFT/blit_texture_complex_to_real_r.comp"));
	std::shared_ptr<ComputeProgram> cp_blit_texture_complex_to_real_i	= std::make_shared<ComputeProgram>(Shader("../CTAnalyzer/Source/GLSL/Compute/FFT/blit_texture_complex_to_real_i.comp"));

	std::shared_ptr<ComputeProgram> cp_conjugate_texture				= std::make_shared<ComputeProgram>(Shader("../CTAnalyzer/Source/GLSL/Compute/FFT/conjugate_texture.comp"));
	std::shared_ptr<ComputeProgram> cp_multiply_complex_texture			= std::make_shared<ComputeProgram>(Shader("../CTAnalyzer/Source/GLSL/Compute/FFT/multiply_complex_texture.comp"));
	std::shared_ptr<ComputeProgram> cp_divide_complex_texture			= std::make_shared<ComputeProgram>(Shader("../CTAnalyzer/Source/GLSL/Compute/FFT/divide_complex_texture.comp"));

	int floor_log2(int n);
	int ceil_log2(int n);

private:
};