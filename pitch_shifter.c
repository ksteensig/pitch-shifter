#define MEOW_FFT_IMPLEMENTATION

#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include "meow_fft.h"

#define WINDOW_SIZE 1024
#define HOP_SIZE 256

char audio_file_name[256];
FILE *audio_file;
float *audio_data;
uint32_t audio_data_len;

float step_size;
float alpha;

uint32_t hop_out;
float hop_buf[HOP_SIZE * 3];

float hann_window[WINDOW_SIZE] = {0};

float ring_buffer[WINDOW_SIZE + HOP_SIZE];
float fft_buffer_in[WINDOW_SIZE];
Meow_FFT_Complex fft_buffer_out[WINDOW_SIZE];
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
  alpha = pow(2, step_size / 12);
  hop_out = round(HOP_SIZE * alpha);

  for (int i = 0; i < WINDOW_SIZE; i++) {
    hann_window[i] = 0.5 * (1 - cos(2 * M_PI * i / (WINDOW_SIZE - 1)));
  }

  workset_bytes = meow_fft_generate_workset_real(WINDOW_SIZE, NULL);
  fft_real = (Meow_FFT_Workset_Real *)malloc(workset_bytes);
  meow_fft_generate_workset_real(WINDOW_SIZE, fft_real);
}

static inline float meow_abs(float r, float j) {
  return abs(sqrt(r * r + j * j));
}

static inline float meow_angle(float r, float j) { return atan(r / j); }

// bit reversal from the FFT is not important as long as
// coefficients are bit reversed too
void pitch_shift() {
  /* Analysis */

  for (int i = 0; i < WINDOW_SIZE; i++) {
    fft_buffer_in[i] *= hann_window[i] * 0.7071;
  }

  // simple fft
  meow_fft_real(fft_real, fft_buffer_in, fft_buffer_out);

  float r, j;

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
    delta_phi[i] = fmod(delta_phi[i] + M_PI, 2 * M_PI) - M_PI;

    prev_phase_frame[i] = phase_frame[i];

    true_freq[i] = (delta_phi_const * i) + delta_phi[i] * (1.0 / HOP_SIZE);

    cumulative_phase[i] += hop_out * true_freq[i];

    /* Synthesis */
  }
}

int main(int argc, char *argv[]) {
  configure();

  pitch_shift();

  /*
  if (argc >= 3) {
    sprintf(argv[1], "%s", audio_file_name);
    sprintf(argv[2], "%lf", step_size);

    printf("%s, %lf\n", audio_file_name, step_size);
  } else {
    printf("Not enough arguments\n");
    printf(
        "First argument is the audio file, second argument is the step size\n");
    return -1;
  }*/

  free(fft_real);

  return 0;
}