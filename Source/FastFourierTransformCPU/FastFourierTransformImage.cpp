#include <vector>
#include <functional>
#include <iostream>

#include "Debuger.h"
#include "FastFourierTransformCPU/FastFourierTransform.h"
#include "FastFourierTransformCPU/ParallelFastFourierTransform.h"
#include "FastFourierTransformCPU/FastFourierTransformImage.h"

void FFT::image_fft_1d_radix2(Image& image) {
	unsigned char* image_data = image.get_image_data();
	std::vector<PFFT::complex> fft_signal(image.get_width());

	ASSERT(image.get_channel_count() == 1);

	int width = image.get_width();
	int height = image.get_height();
	int stride = image.get_byte_per_channel();

	for (int y = 0; y < image.get_height(); y++) {

		fft_signal.clear();
		for (int x = 0; x < image.get_width(); x++) {
			if (image.get_byte_per_channel() == 1)
				fft_signal.push_back(((unsigned char*)image_data)[y * width + x]);
			else if (image.get_byte_per_channel() == 2)
				fft_signal.push_back(((unsigned short*)image_data)[y * width + x]);
			else if (image.get_byte_per_channel() == 4)
				fft_signal.push_back(((unsigned int*)image_data)[y * width + x]);
			else {
				ASSERT(false);
			}
		}

		std::function<double(double)> ram_lak_filter =		[](double x) { return std::abs(x); };
		std::function<double(double)> shepp_logan_filter =	[](double x) { return std::sin(std::abs(x)); };
		std::function<double(double)> cosine_filter =		[](double x) { return std::abs(sin(PI * x)) / PI; };

		PFFT::fft_radix2(fft_signal);
		PFFT::fft_shift(fft_signal);
		
		int size = fft_signal.size();
		
		for (int i = 0; i < fft_signal.size(); i++) {
			double filter_value = shepp_logan_filter((float)i / size - 0.5f) * size;
			fft_signal[i] = PFFT::mult(filter_value, fft_signal[i]);
		}
		
		//for (int i = 0; i < size / 2; i++) {
		//	fft_signal[i] = PFFT::mult(cosine_filter((i) / (float)size) * size, fft_signal[i]);
		//}
		//for (int i = size / 2; i < size; i++) {
		//	fft_signal[i] = PFFT::mult(cosine_filter((size - i) / (float)size) * size, fft_signal[i]);
		//}
		PFFT::inverse_fft_shift(fft_signal);
		
		// Length variables
		//size_t n = fft_signal.size();
		//int log2_size = PFFT::_floor_log2(n);
		//
		//if (static_cast<size_t>(1U) << log2_size != n)
		//	throw std::domain_error("Length is not a power of 2");
		//
		//std::vector<PFFT::complex> exp_table(n / 2);
		//size_t i;
		//for (i = 0; i < n / 2; i++)
		//	exp_table[i] = PFFT::polar(1.0f, -2 * PFFT::PI * i / n);
		//
		//for (i = 0; i < n; i++) {
		//	size_t j = PFFT::_reverse_bits(i, log2_size);
		//	if (j > i)
		//		std::swap(fft_signal[i], fft_signal[j]);
		//}

		//PFFT::_parallel_fft_reverse_bit_order(fft_signal);

		double max = -(1 << 16);
		for (int i = 0; i < fft_signal.size(); i++) {
			max = std::max(fft_signal[i].r , max);
		}
		std::cout << max << std::endl;

		for (int x = 0; x < image.get_width(); x++) {
			if (image.get_byte_per_channel() == 1)
				((unsigned char*)image_data)[y * width + x] = std::max(std::min(fft_signal[x].r / 2 + 128 - 1, 255.0), 0.0);
			else if (image.get_byte_per_channel() == 2)
				((unsigned short*)image_data)[y * width + x] = std::max(std::min(fft_signal[x].r / 2 + 128 * 256 - 1, 256.0 * 256.0 - 1), 0.0);
			else if (image.get_byte_per_channel() == 4)
				((unsigned int*)image_data)[y * width + x] = std::max(std::min(fft_signal[x].r / 2 - 128 * 256 * 256 * 256 - 1, 256.0 * 256.0 * 256.0 * 256.0 - 1), 0.0);
			else {
				ASSERT(false);
			}
		}
	}
}