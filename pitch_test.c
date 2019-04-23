#define MEOW_FFT_IMPLEMENTATION

#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include "meow_fft.h"

#define WINDOW_SIZE 1024
#define HOP_SIZE 256
#define STEP_SIZE 1

char audio_file_name[] = "input_audio.txt";
FILE *audio_file;
float *audio_data;
uint32_t audio_data_len;
int eof = 0;

FILE *output_file;

float step_size;
float alpha;

uint32_t hop_out;
float hop_buf[HOP_SIZE * 3];

uint32_t in_ring_buf_len = WINDOW_SIZE + HOP_SIZE;
uint32_t in_ring_buf_ptr = 0;
float input_ring_buf[WINDOW_SIZE + HOP_SIZE];

// uint32_t out_ring_buf_len = 2300;
uint32_t out_ring_buf_len = 1280;
uint32_t out_ring_buf_ptr = 0;
float output_ring_buf[1800];

const float window_scale_factor = 0.7071;
float hann_window[WINDOW_SIZE] = {0};

float ring_buffer[WINDOW_SIZE + HOP_SIZE];
float fft_buffer_in[WINDOW_SIZE];
Meow_FFT_Complex fft_buffer_out[WINDOW_SIZE];
Meow_FFT_Complex fft_buffer_temp[WINDOW_SIZE];
size_t workset_bytes;
Meow_FFT_Workset_Real *fft_real;

float magnitude_frame[WINDOW_SIZE];
float phase_frame[WINDOW_SIZE];
float prev_phase_frame[WINDOW_SIZE];
float cumulative_phase[WINDOW_SIZE] = {0};

const float delta_phi_const = (2 * M_PI) / WINDOW_SIZE;
float delta_phi[WINDOW_SIZE];

float true_freq[WINDOW_SIZE];

void configure() {
  audio_file = fopen(audio_file_name, "r");
  output_file = fopen("output_audio.txt", "w");

  alpha = pow(2, STEP_SIZE / 12);
  hop_out = round(HOP_SIZE * alpha);

  for (int i = 0; i < WINDOW_SIZE; i++) {
    hann_window[i] = pow(sin((M_PI * i) * (1.0 / WINDOW_SIZE)), 2);
  }

  workset_bytes = meow_fft_generate_workset_real(WINDOW_SIZE, NULL);
  fft_real = (Meow_FFT_Workset_Real *)malloc(workset_bytes);
  meow_fft_generate_workset_real(WINDOW_SIZE, fft_real);
}

static inline float meow_abs(float r, float j) { return sqrt(r * r + j * j); }

static inline float meow_angle(float r, float j) {
  double x = atan2(j, r);
  // printf("%lf, %lf, %lf\n", x, r, j);
  return x;
}

// bit reversal from the FFT is not important as long as
// coefficients are bit reversed too
void time_stretch() {
  /* Analysis */

  for (int i = 0; i < WINDOW_SIZE; i++) {
    // printf("%lf\n", fft_buffer_in[i]);
    fft_buffer_in[i] *= hann_window[i] * window_scale_factor;
  }

  // simple fft
  meow_fft_real(fft_real, fft_buffer_in, fft_buffer_out);

  float r, j;
  float expr, expj;

  // first part of for-loop is still in the analysis part
  // this is done to save WINDOW_SIZE amount of jumps
  for (int i = 0; i < WINDOW_SIZE; i++) {
    r = fft_buffer_out[i].r;
    j = fft_buffer_out[i].j;

    // abs mnemonic
    magnitude_frame[i] = meow_abs(r, j);

    // some atan intrinsic on the TMS320C55X can do this
    phase_frame[i] = meow_angle(r, j);

    /* Processing */
    delta_phi[i] =
        phase_frame[i] - prev_phase_frame[i] - (HOP_SIZE * delta_phi_const * i);

    prev_phase_frame[i] = phase_frame[i];

    // page 87 in the following pdf shows fmod exists for the c55x in assembly
    // https://www.eecs.umich.edu/courses/eecs452/Labs/Docs/C55x_Assmbly_Lang_guide.pdf
    delta_phi[i] = fmod(delta_phi[i] + M_PI, 2 * M_PI) - M_PI;

    true_freq[i] = (delta_phi_const * i) + delta_phi[i] * (1.0 / HOP_SIZE);

    cumulative_phase[i] += hop_out * true_freq[i];

    /* Synthesis */

    expr = cos(cumulative_phase[i]);
    expj = sin(cumulative_phase[i]);

    fft_buffer_out[i].r = magnitude_frame[i] * expr;
    fft_buffer_out[i].j = magnitude_frame[i] * expj;
  }

  // old output buffer is now the input and vice versa
  meow_fft_real_i(fft_real, fft_buffer_out, fft_buffer_temp, fft_buffer_in);

  for (int i = 0; i < WINDOW_SIZE; i++) {
    fft_buffer_in[i] = fft_buffer_in[i] * hann_window[i] * window_scale_factor *
                       (1.0 / WINDOW_SIZE);
  }
}

int main(int argc, char *argv[]) {
  configure();

  float p;

  for (int i = 0; i < 6000; i++) {
    fscanf(audio_file, "%f\n", &p);
  }

  for (int i = 0; i < WINDOW_SIZE; i++) {
    fscanf(audio_file, "%f\n", &fft_buffer_in[i]);
  }

  time_stretch();

  for (int i = 0; i < WINDOW_SIZE; i++) {
    fprintf(output_file, "%f\n", fft_buffer_in[i]);
  }

  fclose(audio_file);
  fclose(output_file);

  free(fft_real);

  return 0;
}