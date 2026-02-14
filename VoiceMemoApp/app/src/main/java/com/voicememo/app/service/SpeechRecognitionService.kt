package com.voicememo.app.service

import android.content.Context
import android.content.Intent
import android.os.Bundle
import android.speech.RecognitionListener
import android.speech.RecognizerIntent
import android.speech.SpeechRecognizer
import kotlinx.coroutines.flow.MutableStateFlow
import kotlinx.coroutines.flow.StateFlow
import kotlinx.coroutines.flow.asStateFlow

class SpeechRecognitionService(private val context: Context) {

    sealed class RecognitionState {
        data object Idle : RecognitionState()
        data object Listening : RecognitionState()
        data class Partial(val text: String) : RecognitionState()
        data class Result(val text: String) : RecognitionState()
        data class Error(val message: String) : RecognitionState()
    }

    private val _state = MutableStateFlow<RecognitionState>(RecognitionState.Idle)
    val state: StateFlow<RecognitionState> = _state.asStateFlow()

    private var speechRecognizer: SpeechRecognizer? = null
    private val accumulatedText = StringBuilder()

    private val recognitionIntent = Intent(RecognizerIntent.ACTION_RECOGNIZE_SPEECH).apply {
        putExtra(RecognizerIntent.EXTRA_LANGUAGE_MODEL, RecognizerIntent.LANGUAGE_MODEL_FREE_FORM)
        putExtra(RecognizerIntent.EXTRA_LANGUAGE, "nl-NL")
        putExtra(RecognizerIntent.EXTRA_PARTIAL_RESULTS, true)
        putExtra(RecognizerIntent.EXTRA_MAX_RESULTS, 1)
        // Prefer offline recognition
        putExtra(RecognizerIntent.EXTRA_PREFER_OFFLINE, true)
    }

    fun startListening() {
        if (speechRecognizer == null) {
            speechRecognizer = SpeechRecognizer.createSpeechRecognizer(context)
            speechRecognizer?.setRecognitionListener(createListener())
        }
        _state.value = RecognitionState.Listening
        speechRecognizer?.startListening(recognitionIntent)
    }

    fun stopListening() {
        speechRecognizer?.stopListening()
    }

    fun getAccumulatedText(): String = accumulatedText.toString()

    fun clearAccumulatedText() {
        accumulatedText.clear()
    }

    fun destroy() {
        speechRecognizer?.destroy()
        speechRecognizer = null
    }

    private fun createListener() = object : RecognitionListener {
        override fun onReadyForSpeech(params: Bundle?) {
            _state.value = RecognitionState.Listening
        }

        override fun onBeginningOfSpeech() {}

        override fun onRmsChanged(rmsdB: Float) {}

        override fun onBufferReceived(buffer: ByteArray?) {}

        override fun onEndOfSpeech() {}

        override fun onError(error: Int) {
            val message = when (error) {
                SpeechRecognizer.ERROR_AUDIO -> "Audiofout"
                SpeechRecognizer.ERROR_CLIENT -> "Clientfout"
                SpeechRecognizer.ERROR_INSUFFICIENT_PERMISSIONS -> "Geen microfoontoestemming"
                SpeechRecognizer.ERROR_NETWORK -> "Netwerkfout"
                SpeechRecognizer.ERROR_NETWORK_TIMEOUT -> "Netwerk timeout"
                SpeechRecognizer.ERROR_NO_MATCH -> "Geen spraak herkend"
                SpeechRecognizer.ERROR_RECOGNIZER_BUSY -> "Herkenner bezet"
                SpeechRecognizer.ERROR_SERVER -> "Serverfout"
                SpeechRecognizer.ERROR_SPEECH_TIMEOUT -> "Geen spraak gedetecteerd"
                else -> "Onbekende fout ($error)"
            }
            // For NO_MATCH and SPEECH_TIMEOUT, auto-restart to enable continuous recording
            if (error == SpeechRecognizer.ERROR_NO_MATCH || error == SpeechRecognizer.ERROR_SPEECH_TIMEOUT) {
                speechRecognizer?.startListening(recognitionIntent)
                return
            }
            _state.value = RecognitionState.Error(message)
        }

        override fun onResults(results: Bundle?) {
            val matches = results?.getStringArrayList(SpeechRecognizer.RESULTS_RECOGNITION)
            val text = matches?.firstOrNull() ?: ""
            if (text.isNotBlank()) {
                if (accumulatedText.isNotEmpty()) {
                    accumulatedText.append(" ")
                }
                accumulatedText.append(text)
            }
            _state.value = RecognitionState.Result(accumulatedText.toString())
            // Auto-restart for continuous listening
            speechRecognizer?.startListening(recognitionIntent)
        }

        override fun onPartialResults(partialResults: Bundle?) {
            val matches = partialResults?.getStringArrayList(SpeechRecognizer.RESULTS_RECOGNITION)
            val partial = matches?.firstOrNull() ?: ""
            if (partial.isNotBlank()) {
                val fullText = if (accumulatedText.isNotEmpty()) {
                    "${accumulatedText} $partial"
                } else {
                    partial
                }
                _state.value = RecognitionState.Partial(fullText)
            }
        }

        override fun onEvent(eventType: Int, params: Bundle?) {}
    }
}
