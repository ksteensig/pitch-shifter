#define MEOW_FFT_IMPLEMENTATION

#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include "meow_fft.h"

#define WINDOW_SIZE 1024
#define HOP_SIZE 256
#define STEP_SIZE 2

char audio_file_name[] = "input_audio.txt";
FILE *audio_file;
float *audio_data;
uint32_t audio_data_len;
int eof = 0;

FILE *output_file;

float alpha;
uint32_t hop_out;

uint32_t in_ring_buf_len = WINDOW_SIZE + HOP_SIZE;
uint32_t in_ring_buf_ptr = 0;
float input_ring_buf[WINDOW_SIZE + HOP_SIZE];

uint32_t out_ring_buf_len = 2300;
uint32_t out_ring_buf_ptr = 0;
float output_ring_buf[2300];

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

float floatMod(float a, float b)
{
    return (a - b * floor(a * 1/b));
}

void configure() {
  audio_file = fopen(audio_file_name, "r");
  output_file = fopen("output_audio.txt", "w");

  alpha = pow(2, STEP_SIZE / 12.0);
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
  for (int i = 0; i < WINDOW_SIZE; i++) {
    // printf("%f\n", input_ring_buf[in_ring_buf_ptr]);
    fft_buffer_in[i] = input_ring_buf[in_ring_buf_ptr];
    in_ring_buf_ptr = ++in_ring_buf_ptr % in_ring_buf_len;
  }

  /* Analysis */

  for (int i = 0; i < WINDOW_SIZE; i++) {
    // printf("%lf\n", fft_buffer_in[i]);
    fft_buffer_in[i] = fft_buffer_in[i] * hann_window[i] * window_scale_factor;
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

    // page 87 in the following pdf shows fmod exists for the c55x in assembly
    // https://www.eecs.umich.edu/courses/eecs452/Labs/Docs/C55x_Assmbly_Lang_guide.pdf
    delta_phi[i] = floatMod(delta_phi[i] + M_PI, 2*M_PI) - M_PI;

    prev_phase_frame[i] = phase_frame[i];

    true_freq[i] = (delta_phi_const * i) + delta_phi[i] * (1.0 / HOP_SIZE);

    cumulative_phase[i] += hop_out * true_freq[i];
    cumulative_phase[i] = floatMod(cumulative_phase[i], 2*M_PI);

    /* Synthesis */

    expr = cos(cumulative_phase[i]);
    expj = sin(cumulative_phase[i]);

    fft_buffer_out[i].r = magnitude_frame[i] * expr;
    fft_buffer_out[i].j = magnitude_frame[i] * expj;
  }

  // old output buffer is now the input and vice versa
  meow_fft_real_i(fft_real, fft_buffer_out, fft_buffer_temp, fft_buffer_in);

  uint32_t start_out = out_ring_buf_ptr;

  for (int i = 0; i < WINDOW_SIZE; i++) {
    fft_buffer_in[i] = fft_buffer_in[i] * hann_window[i] * window_scale_factor *
                       (1.0 / WINDOW_SIZE);
    output_ring_buf[out_ring_buf_ptr] += fft_buffer_in[i];
    out_ring_buf_ptr = ++out_ring_buf_ptr % out_ring_buf_len;
  }

  // set old zone to zero
  for (int i = 0; i < hop_out; i++) {
    output_ring_buf[(out_ring_buf_ptr + i) % out_ring_buf_len] = 0;
  }

  // move one synthetic hop size in the output ring buffer
  out_ring_buf_ptr = (start_out + hop_out) % out_ring_buf_len;
}

void resample() {
  uint32_t x1, x2;
  float a;

  for (int i = 0; i < WINDOW_SIZE; i++) {
    // linear interpolation
    x1 = (uint32_t)(out_ring_buf_ptr + alpha * i) % out_ring_buf_len;
    x2 = (x1 + 1) % out_ring_buf_len;
    a = output_ring_buf[x2] - output_ring_buf[x1];
    fft_buffer_in[i] = output_ring_buf[x1] + alpha * a;
  }
}

inline void get_sample() {
  if (fscanf(audio_file, "%f\n", &input_ring_buf[in_ring_buf_ptr]) == EOF) {
    eof = 1;
    input_ring_buf[in_ring_buf_ptr] = 0;
  }

  // printf("%f\n", input_ring_buf[in_ring_buf_ptr]);

  in_ring_buf_ptr = ++in_ring_buf_ptr % in_ring_buf_len;
}

inline void set_samples() {
  for (int i = 0; i < WINDOW_SIZE; i++) {
    fprintf(output_file, "%f\n", fft_buffer_in[i]);
  }
}

int main(int argc, char *argv[]) {
  configure();

  uint32_t start_in;
  uint32_t start_out;

  while (!eof) {
    start_out = out_ring_buf_ptr;
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < HOP_SIZE; j++) {
        get_sample();
      }

      start_in = in_ring_buf_ptr;
      in_ring_buf_ptr = (in_ring_buf_ptr + HOP_SIZE) % in_ring_buf_len;

      time_stretch();

      in_ring_buf_ptr = start_in;
    }

    out_ring_buf_ptr = start_out;

    resample();

    out_ring_buf_ptr =
        (start_out + (uint32_t)round(4 * hop_out)) % out_ring_buf_len;

    set_samples();
  }

  fclose(audio_file);
  fclose(output_file);

  free(fft_real);

  return 0;
}