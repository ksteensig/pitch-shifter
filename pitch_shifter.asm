; Note that this whole ASM file is currently pseudo code

TWO_PI_Q5_11        .set 0011001001000100b
PI_Q5_11            .set 0001100100100010b

WINDOW_SIZE         .set 1024
INV_WINDOW_SIZE_Q15 .set 0000000000100000b
HOP_SIZE		        .set 256
INV_HOP_SIZE_Q15    .set 0000000010000000b
STEP_SIZE           .set 1
WINDOW_SCALE_FACTOR .set 0x59BA ;Q15 0.701    ; Constant that comes from sqrt(((WINDOW_SIZE/HOP_SIZE)/2))
DELTA_PHI_CONST     .set 0x00C9 ; Q15 0.00613  ; approx 2pi/WINDOW_SIZE

HANN_WINDOW_START       .set 0x0000 ; TBD
FFT_BUF_START           .set 0x0000 ; TBD
MAGNITUDE_FRAME_START   .set 0x0000 ; TBD
PHASE_FRAME_START       .set 0x0000 ; TBD
PREV_PHASE_FRAME_START  .set 0x0000 ; TBD
DELTA_PHI_START         .set 0x0000 ; TBD
BIT_REV_INDEX_START     .set 0x0000 ; TBD
TRUE_FREQ_START         .set 0x0000 ; TBD
CUMMULATIVE_PHASE_START .set 0x0000 ; TBD

HOP_OUT                 .set 0x0000 ; TBD
ALPHA                   .set 0x0000 ; TBD

ATAN2_R                 .set 0x0000 ; TBD
ATAN2_J                 .set 0x0000 ; TBD

FFT_BUF         .usect "PITCH_SHIFT_FFT_BUF",      2048
INPUT_RING_BUF  .usect "PITCH_SHIFT_RING_BUF_IN",  1200 ; 5*HOP_SIZE
OUTPUT_RING_BUF .usect "PITCH_SHIFT_RING_BUF_OUT", 2300 ; 7*HOP_SIZE*2^(4/12) rounded up to nearest 100

HANN_WINDOW       .usect "HANN_WINDOW",             1024
MAGNITUDE_FRAME   .usect "MAGNITUDE_FRAME",         1024
PHASE_FRAME       .usect "PHASE_FRAME",             1024
PREV_PHASE_FRAME  .usect "PREV_PHASE_FRAME",        1024
BIT_REV_INDEX     .usect "BIT_REVERSED_INDEX",      1024
DELTA_PHI         .usect "DELTA_PHI",               1024
TRUE_FREQ         .usect "TRUE_FREQ",               1024
CUMMULATIVE_PHASE .usect "CUMMULATIVE_PHASE",       1024


_cifft_SCALE:
  RET

_cfft_SCALE:
  RET

_sqrt_16:
  RET

_atan16:
  RET

** T0 = input (Q5.11)
** T1 = output (Q5.11)
_mod_2pi: .macro
  NOP
  .endm
  

_time_stretch:
  MOV #HANN_WINDOW_START, AR0           ; Need a pointer to the hanning window
  MOV #1023, BRC0                       ; RTPBLOCAL will loop BRC0+1 times
  RPTBLOCAL ps_loop1_end                ; Loop with WINDOW_SIZE iterations to move new window into work buffer
  MOV *AR4+, AC0                        ; Assume AR4 is the input FFT_BUF
  MOV *AR0, AC1     
  MPY AC1, AC0                          ; Window the sample
  MPY #WINDOW_SCALE_FACTOR, AC0
  ADD #1, AR0
ps_loop1_end: ADD #1, AR4
                                        ; FFT the frame
  SUB #2048, AR4
  MOV #1024, T0
  CALL _cfft_SCALE                      ; Best performance with twiddle in ROM and FFT_BUF in DARAM
                                        
                                        ; Calculate magnitude frame
  MOV #MAGNITUDE_FRAME_START, AR0
  MOV #1023, BRC0
  RPTBLOCAL mag_frame_loop_end
  MOV *AR4+, AC0
  MOV *AR4+, AC1
  SQR AC0
  SQR AC1
  ADD AC1, AC0
mag_frame_loop_end: MOV AC0, *AR0+
  CALL _sqrt_16                         ; Call sqrt on all magnitude frame members

                                        ; Calculate the phase frame, the angles for each sample
  MOV #ATAN2_R, AR0
  MOV #ATAN2_J, AR1
  SUB #2048, AR4
  MOV #1023, BRC0                       ; Another 1024 iterations
  RPTBLOCAL ps_loop2_end
  MOV *AR4+, *AR0+
ps_loop2_end: MOV *AR4+, *AR1+
  MOV #ATAN2_R, AR0                     ; Set pointers back to beginning
  MOV #ATAN2_J, AR1
  MOV #PHASE_FRAME_START, AR2
  MOV #1024, T0
  CALL _atan16
                                        ; Calculate delta phi
  MOV 1023, BRC0
  MOV #DELTA_PHI_START, AR0
  MOV #PHASE_FRAME_START, AR1
  MOV #PREV_PHASE_FRAME_START, AR2
  MOV #BIT_REV_INDEX, AR3

  RPTBLOCAL ps_loop3_end
  MOV *AR0, AC0
  MOV #HOP_SIZE, AC0
  MPY #DELTA_PHI_CONST, AC0
  MPY *AR3+, AC0
  NEG AC0
  ADD *AR1, AC0
  SUB *AR2, AC0
  ADD #PI_Q5_11, AC0
  _mod_2pi
  SUB #PI_Q5_11, AC0
  MOV *AR1+, *AR2+
  MOV AC0, T0
ps_loop3_end: MOV T1, *AR0+

  MOV #1023, BRC0
  SUB 1024, AR0                         ; Go back to beginning of delta phi vector
  MOV #TRUE_FREQ_START, AR1
  SUB #1024, AR3                        ; Go back to beginning of bit reversed index vector
  
  RPTBLOCAL ps_loop4_end
  MOV #DELTA_PHI_CONST, AC0
  MPY *AR3+, AC0
  MOV *AR0+, AC1
  MPY #INV_HOP_SIZE_Q15, AC1
  ADD AC1, AC0
ps_loop4_end: MOV AC1, *AR1+

                                        ; Calculate the cumulative phase
  MOV #1023, BRC0
  MOV #CUMMULATIVE_PHASE_START, AR0
  SUB #1024, AR1                        ; Go back to beginning of true freq vector
  MOV #HOP_OUT, AR2
  RPTBLOCAL ps_loop5_end
  MOV *AR2, AC0
  MPY *AR1+, AC0
  MOV AC0, T0
  _mod_2pi
ps_loop5_end: MOV T1, *AR0+
                                        ; Phase shift the FFT'd buffer, assume AR7 pointsa at sine wave
  MOV #1023, BRC0
  SUB #1024, AR0                        ; Go back to beginning of cummulative phase
  SUB #2048, AR4                        ; Go back to beginning of the FFT buffer

  RPTBLOCAL ps_loop6_end
  ;MOV *AR0+, T0
  ;MOV *AR7(T0), 
ps_loop6_end: NOP

  MOV #FFT_BUF_START, AR0
  MOV #1024, T0
  CALL _cifft_SCALE

  MOV #1023, BRC0
  SUB #2048, AR4
  MOV #HANN_WINDOW_START, AR0

  RPTBLOCAL ps_loop7_end
  MOV *AR4+, AC0
  MPY *AR0+, AC0
  MPY #WINDOW_SCALE_FACTOR, AC0
  MPY #INV_WINDOW_SIZE_Q15, AC0
  ADD *AR7, AC0
ps_loop7_end: MOV AC0, *AR7+
                                        ; Zero old values so it won't keep accumulating
  MOV #HOP_OUT, BRC0
  RPTBLOCAL ps_loop8_end
ps_loop8_end: MOV #0, *AR7+

  SUB #WINDOW_SIZE, AR7                 ; Only subtract the window size, so it will have kipped hop out size
  RET

_resample:
  MOV #HOP_OUT, AC0
  MPY #4, AC0
  SUB AC0, AR7

  MOV #1023, BRC0
  ; Make AR0 point at the output buffer
  RPTBLOCAL ps_loop9_end
  MOV *AR7+, AC1
  MOV *AR7+, AC0
  SUB AC1, AC0
  MPY #ALPHA, AC1
  ADD AC1, AC0
  MOV AC0, *AR0+
ps_loop9_end: SUB #1, AR7               ; Go one back so it will only loop 1024 times over AR7

  MOV HOP_OUT, AC0
  MPY #4, AC0
  ADD AC0, AR7
  RET