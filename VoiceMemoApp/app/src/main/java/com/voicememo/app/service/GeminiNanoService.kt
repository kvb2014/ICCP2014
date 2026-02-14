package com.voicememo.app.service

import android.content.Context
import com.google.ai.edge.aicore.GenerativeModel
import com.google.ai.edge.aicore.generationConfig

/**
 * Service that uses Gemini Nano (on-device via Google AI Edge) to
 * restructure a voice memo transcript into a structured markdown note
 * following a user-provided template and instructions.
 */
class GeminiNanoService(private val context: Context) {

    private var generativeModel: GenerativeModel? = null

    private fun getOrCreateModel(): GenerativeModel {
        return generativeModel ?: GenerativeModel(
            generationConfig {
                temperature = 0.2f
                topK = 16
                maxOutputTokens = 4096
            }
        ).also { generativeModel = it }
    }

    /**
     * Processes a transcript using Gemini Nano, structuring it according to
     * the provided template and instructions.
     *
     * @param transcript The raw voice memo transcript
     * @param templateContent The markdown template to follow
     * @param templateInstructions Additional instructions for the LLM
     * @param existingNote If continuing an existing note, its current content
     * @return The structured markdown note
     */
    suspend fun processTranscript(
        transcript: String,
        templateContent: String,
        templateInstructions: String,
        existingNote: String? = null
    ): String {
        val model = getOrCreateModel()

        val prompt = buildString {
            append("Je bent een assistent die spraaknotities omzet naar gestructureerde markdown aantekeningen.\n\n")

            if (templateInstructions.isNotBlank()) {
                append("## Instructies\n")
                append(templateInstructions)
                append("\n\n")
            }

            append("## Template\n")
            append("Gebruik het volgende markdown-template als structuur voor de aantekening:\n\n")
            append(templateContent)
            append("\n\n")

            if (existingNote != null) {
                append("## Bestaande aantekening\n")
                append("Er is al een eerdere versie van deze aantekening. Vul deze aan met de nieuwe informatie uit het transcript.\n\n")
                append(existingNote)
                append("\n\n")
            }

            append("## Transcript\n")
            append("Hier is het transcript van de spraakmemo:\n\n")
            append(transcript)
            append("\n\n")

            append("## Opdracht\n")
            append("Genereer nu de gestructureerde markdown aantekening op basis van bovenstaand transcript en template. ")
            append("Vul alle velden in het template zo goed mogelijk in op basis van de inhoud van het transcript. ")
            append("Geef alleen de markdown aantekening terug, zonder extra uitleg.")
        }

        val response = model.generateContent(prompt)
        return response.text ?: "Fout: geen resultaat van Gemini Nano"
    }

    fun close() {
        generativeModel?.close()
        generativeModel = null
    }
}
