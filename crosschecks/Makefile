LIGO_FRAME_DURATION=64
LIGO_N_FRAMES=1

INJ_TIME=$(shell ligolw_print SG*xml -t sim_burst -c time_geocent_gps)
START=$(shell python -c "print ${INJ_TIME}-0.5*${LIGO_FRAME_DURATION}")
STOP=$(shell python -c "print ${INJ_TIME}+0.5*${LIGO_FRAME_DURATION}")
DUR=$(shell python -c "print ${STOP}-${START}")

BURSTMDC_LOG=SG361_test_BurstMDC-ER5BURST-Log.txt
INJ_XML=SG361_test_BurstMDC-ER5BURST.xml
INJ_CHANNEL_NAME=FAKE-STRAIN
INJECTION_FRAME_TYPE=SG361_test_BurstMDC_ER5BURST
INJECTION_DEST=./

LIGO_FRAME_DURATION=64
LIGO_N_FRAMES=1

all: H1_frames L1_frames V1_frames

#
# Injection tables
#
timeslide.xml:
	ligolw_tisi -i H1=0:0:1 -i L1=0:0:1 -i V1=0:0:1 > timeslide.xml

SG361_test_BurstMDC-ER5BURST.xml: timeslide.xml
	ligolw_burstmdc_to_simburst \
		-t timeslide.xml \
		${BURSTMDC_LOG} \
		-o ${INJ_XML}

H1_frames: ${INJ_XML}
	gstlal_fake_frames --gps-start-time ${START} --gps-end-time ${STOP} --injections ${INJ_XML} --channel H1=${INJ_CHANNEL_NAME} --verbose --data-source silence --output-channel ${INJ_CHANNEL_NAME} --frame-type H1_${INJECTION_FRAME_TYPE} --output-path ${INJECTION_DEST} --frame-duration ${LIGO_FRAME_DURATION} --frames-per-file ${LIGO_N_FRAMES}

L1_frames: ${INJ_XML}
	gstlal_fake_frames --gps-start-time ${START} --gps-end-time ${STOP} --injections ${INJ_XML} --channel L1=${INJ_CHANNEL_NAME} --verbose --data-source silence --output-channel ${INJ_CHANNEL_NAME} --frame-type L1_${INJECTION_FRAME_TYPE} --output-path ${INJECTION_DEST} --frame-duration ${LIGO_FRAME_DURATION} --frames-per-file ${LIGO_N_FRAMES}

V1_frames: ${INJ_XML}
	gstlal_fake_frames --gps-start-time ${START} --gps-end-time ${STOP} --injections ${INJ_XML} --channel V1=${INJ_CHANNEL_NAME} --verbose --data-source silence --output-channel ${INJ_CHANNEL_NAME} --frame-type V1_${INJECTION_FRAME_TYPE} --output-path ${INJECTION_DEST} --frame-duration ${LIGO_FRAME_DURATION} --frames-per-file ${LIGO_N_FRAMES}

clean:
	rm -rf *H1_SG361_test_BurstMDC_ER5BURST*
	rm -rf *L1_SG361_test_BurstMDC_ER5BURST*
	rm -rf *V1_SG361_test_BurstMDC_ER5BURST*
	rm -rf ${INJ_XML}
	rm -rf timeslide.xml

