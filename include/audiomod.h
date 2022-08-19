/*! \mainpage audiomod - sound modification project
 *
 * \section Introduction
 * This project contains a lot of well designed sound modification tools, or digital audio effects (DAFx), including:
 * \li phasevocoder models (time-stretching, pitch-shifting, vocoder, gender changing)
 * \li delay line models (delay, chorus, flanger, vibrato)
 * \li distortion models
 * \li dynamics models (compressor, limiter, dynamic filter)
 * \li filtering models (autowah, wahwah, phaser)
 * \li modulation models (ring modulation, tremolo)
 * 
 * This project is created and maintained by tangkk
 *
 * \section Workflow
 *
 * Please refer to [Related Pages](pages.html) for more details on usage
 * 
 */

#pragma once

#include "dafx/phasevocoder.h"
#include "dafx/chorus.h"
#include "dafx/delay.h"
#include "dafx/flanger.h"
#include "dafx/vibrato.h"

#include "dafx/distortion.h"

#include "dafx/compressor.h"
#include "dafx/limiter.h"
#include "dafx/dynamicfilter.h"

#include "dafx/autowah.h"
#include "dafx/phaser.h"

#include "dafx/ringmod.h"
#include "dafx/tremolo.h"

#include "dafx/reverb.h"

#include "dafx/loudnessmeter.h"

#include "dafx/biquadfilter.h"
#include "dafx/equalizer.h"

#include "dafx/gain.h"

// #if defined(__APPLE__) or defined(__LINUX__)
// #include "analyzer/vad.h"
#include "analyzer/envelope.h"
// #endif