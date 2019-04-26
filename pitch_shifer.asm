; Note that this whole ASM file is currently pseudo code

WINDOW_SIZE         .set 1024
HOP_SIZE		        .set 256
STEP_SIZE           .set 1
WINDOW_SCALE_FACTOR .set 0.701  ; Constant that comes from sqrt(((WINDOW_SIZE/HOP_SIZE)/2))
DELTA_PHI_CONST     .set (2 * M_PI) / WINDOW_SIZE;

FFT_BUF         .usect "PITCH_SHIFT_FFT_BUF",      1024
IFFT_BUF        .usect "PITCH_SHIFT_IFFT_BUF",     2048
INPUT_RING_BUF  .usect "PITCH_SHIFT_RING_BUF_IN",  1200 ; 5*HOP_SIZE
OUTPUT_RING_BUF .usect "PITCH_SHIFT_RING_BUF_OUT", 2300 ; 7*HOP_SIZE*2^(4/12) rounded up to nearest 100

HANN_WINDOW       .usect "HANN_WINDOW",            1024
MAGNITUDE_FRAME   .usect "MAGNITUDE_FRAME"         1024
PHASE_FRAME       .usect "PHASE_FRAME"             1024
PREV_PHASE_FRAME  .usect "PREV_PHASE_FRAME"        1024
BIT_REV_INDEX     .usect "BIT_REVERSED_INDEX"      1024
DELTA_PHI         .usect "DELTA_PHI"               1024
TRUE_FREQ         .usect "TRUE_FREQ"               1024
CUMMULATIVE_PHASE .usect "CUMMULATIVE_PHASE"       1024

_pitch_shift:
  MOV HANN_WINDOW, T0                   ; Need a pointer to the hanning window
  MOV #1023, BRC0                       ; RTPBLOCAL will loop BRC0+1 times
  RPTBLOCAL ps_loop1_end                ; Loop with WINDOW_SIZE iterations to move new window into work buffer
  MOV *AR4+, AC0                        ; Assume AR4 is the input ringer buffer for the pitch shifter
  MPY T0+, AC0                          ; Window the sample
  MPY WINDOW_SCALE_FACTOR, AC0
ps_loop1_end:
  fft_16 FFT_BUF, IFFT_BUF, WINDOW_SIZE  ; Probably assumes DARAM and isn't called like this

  MOV #0, T0                            ; Create index counter for loop
  MOV #1023 BRC0                        ; Another 1024 iterations
  RTPBLOCAL ps_loop2_end
                                        ; Following instructions calculates abs() of the complex value in the IFFT_BUF
  MOV IFFT_BUF(T0), AC0                 ; AC0 = r, where T0_i = r and T0_i+1 = j from complex numbers
  SQR AC0                               ; AC0 = r^2
  MOV IFFT_BUF(T0+1), AC1               ; AC1 = j
  SQR AC1                               ; AC1 = j^2
  ADD AC1, AC0                          ; AC0 = r^2 + j^2
  sqrt_16 AC0, AC0                      ; sqrt_16 assumes DARAM and isn't called like this
  atan2_16 r, j, PHASE_FRAME, 1         ; calculate the angle with r and j mapped between -pi and pi, also assumes DARAM

                                        ; Calculate phase_frame[i] - prev_phase_frame[i] - (HOP_SIZE * delta_phi_const * i)
  MOV -#HOP_SIZE, AC0
  MPY DELTA_PHI_CONST, AC0
  MPY BIT_REV_INDEX(T0), AC0
  ADD PHASE_FRAME(T0), AC0
  SUB PREV_PHASE_FRAME(T0), AC0
                                        ; Calculate fmod(delta_phi[i] + pi, 2*pi) - pi
  ADD M_PI, AC0
  fmod AC0, 2xM_PI
  SUB M_PI, AC0

  
ps_loop2_end: