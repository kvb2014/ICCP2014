package com.voicememo.app.ui

import android.app.Application
import android.net.Uri
import androidx.lifecycle.AndroidViewModel
import androidx.lifecycle.viewModelScope
import com.voicememo.app.VoiceMemoApplication
import com.voicememo.app.data.model.Memo
import com.voicememo.app.data.model.Template
import com.voicememo.app.service.GeminiNanoService
import com.voicememo.app.service.ObsidianExportService
import com.voicememo.app.service.SpeechRecognitionService
import com.voicememo.app.service.SpeechRecognitionService.RecognitionState
import kotlinx.coroutines.flow.*
import kotlinx.coroutines.launch
import java.text.SimpleDateFormat
import java.util.*

class MainViewModel(application: Application) : AndroidViewModel(application) {

    private val app = application as VoiceMemoApplication
    private val memoRepository = app.memoRepository
    private val templateRepository = app.templateRepository

    val speechService = SpeechRecognitionService(application)
    private val geminiService = GeminiNanoService(application)
    val obsidianService = ObsidianExportService()

    // State
    val memos = memoRepository.allMemos.stateIn(viewModelScope, SharingStarted.Lazily, emptyList())
    val templates = templateRepository.allTemplates.stateIn(viewModelScope, SharingStarted.Lazily, emptyList())

    private val _selectedTemplate = MutableStateFlow<Template?>(null)
    val selectedTemplate: StateFlow<Template?> = _selectedTemplate.asStateFlow()

    private val _selectedVault = MutableStateFlow<String?>(null)
    val selectedVault: StateFlow<String?> = _selectedVault.asStateFlow()

    private val _currentMemoId = MutableStateFlow<Long?>(null)
    val currentMemoId: StateFlow<Long?> = _currentMemoId.asStateFlow()

    private val _currentTranscript = MutableStateFlow("")
    val currentTranscript: StateFlow<String> = _currentTranscript.asStateFlow()

    private val _isProcessing = MutableStateFlow(false)
    val isProcessing: StateFlow<Boolean> = _isProcessing.asStateFlow()

    private val _detailMemo = MutableStateFlow<Memo?>(null)
    val detailMemo: StateFlow<Memo?> = _detailMemo.asStateFlow()

    init {
        // Detect vault on startup
        val vaults = obsidianService.listVaults()
        if (vaults.size == 1) {
            _selectedVault.value = vaults.first()
        }

        // Track speech recognition results
        viewModelScope.launch {
            speechService.state.collect { state ->
                when (state) {
                    is RecognitionState.Result -> {
                        _currentTranscript.value = state.text
                    }
                    is RecognitionState.Partial -> {
                        _currentTranscript.value = state.text
                    }
                    else -> {}
                }
            }
        }
    }

    fun selectTemplate(template: Template) {
        _selectedTemplate.value = template
    }

    fun selectVault(vault: String) {
        _selectedVault.value = vault
    }

    fun startRecording() {
        speechService.startListening()
    }

    fun stopRecording() {
        speechService.stopListening()
        _currentTranscript.value = speechService.getAccumulatedText()
    }

    /**
     * Start recording as a continuation of an existing memo.
     */
    fun continueRecording(memo: Memo) {
        _currentMemoId.value = memo.id
        _currentTranscript.value = memo.transcript
        speechService.clearAccumulatedText()
        // Pre-populate with existing text so new speech appends
        // We'll append manually after recording stops
    }

    fun saveAndProcess() {
        val transcript = speechService.getAccumulatedText().ifBlank { _currentTranscript.value }
        if (transcript.isBlank()) return

        viewModelScope.launch {
            val existingId = _currentMemoId.value
            val dateFormat = SimpleDateFormat("dd-MM-yyyy HH:mm", Locale("nl"))
            val now = System.currentTimeMillis()

            if (existingId != null) {
                // Appending to existing memo
                val existing = memoRepository.getMemoById(existingId)
                if (existing != null) {
                    val updatedTranscript = existing.transcript + "\n\n" + transcript
                    val updatedMemo = existing.copy(
                        transcript = updatedTranscript,
                        updatedAt = now
                    )
                    memoRepository.update(updatedMemo)
                    processWithGemini(updatedMemo)
                }
            } else {
                // New memo
                val title = "Memo ${dateFormat.format(Date(now))}"
                val memo = Memo(
                    title = title,
                    transcript = transcript,
                    templateId = _selectedTemplate.value?.id,
                    createdAt = now,
                    updatedAt = now
                )
                val id = memoRepository.insert(memo)
                val saved = memo.copy(id = id)
                processWithGemini(saved)
            }

            // Reset
            speechService.clearAccumulatedText()
            _currentTranscript.value = ""
            _currentMemoId.value = null
        }
    }

    private fun processWithGemini(memo: Memo) {
        val template = _selectedTemplate.value ?: return
        val vault = _selectedVault.value ?: return

        viewModelScope.launch {
            _isProcessing.value = true
            try {
                val existingNote = memo.processedMarkdown
                val result = geminiService.processTranscript(
                    transcript = memo.transcript,
                    templateContent = template.content,
                    templateInstructions = template.instructions,
                    existingNote = existingNote
                )

                val updatedMemo = memo.copy(
                    processedMarkdown = result,
                    updatedAt = System.currentTimeMillis()
                )
                memoRepository.update(updatedMemo)

                // Export to Obsidian
                val path = obsidianService.saveNote(vault, memo.title, result)
                memoRepository.update(updatedMemo.copy(obsidianPath = path))

                _detailMemo.value = updatedMemo.copy(
                    processedMarkdown = result,
                    obsidianPath = path
                )
            } catch (e: Exception) {
                // If Gemini Nano fails, save transcript as-is in template
                val fallback = "# ${memo.title}\n\n## Transcript\n\n${memo.transcript}"
                val updatedMemo = memo.copy(
                    processedMarkdown = fallback,
                    updatedAt = System.currentTimeMillis()
                )
                memoRepository.update(updatedMemo)
                val path = obsidianService.saveNote(vault, memo.title, fallback)
                memoRepository.update(updatedMemo.copy(obsidianPath = path))
                _detailMemo.value = updatedMemo.copy(
                    processedMarkdown = fallback,
                    obsidianPath = path
                )
            } finally {
                _isProcessing.value = false
            }
        }
    }

    fun loadMemoDetail(memoId: Long) {
        viewModelScope.launch {
            _detailMemo.value = memoRepository.getMemoById(memoId)
        }
    }

    fun reprocessMemo(memoId: Long) {
        viewModelScope.launch {
            val memo = memoRepository.getMemoById(memoId) ?: return@launch
            processWithGemini(memo)
        }
    }

    fun updateMemoMarkdown(memoId: Long, markdown: String) {
        viewModelScope.launch {
            val memo = memoRepository.getMemoById(memoId) ?: return@launch
            val updated = memo.copy(
                processedMarkdown = markdown,
                updatedAt = System.currentTimeMillis()
            )
            memoRepository.update(updated)
            _detailMemo.value = updated

            // Also update in Obsidian
            val vault = _selectedVault.value ?: return@launch
            obsidianService.saveNote(vault, memo.title, markdown)
        }
    }

    fun deleteMemo(memo: Memo) {
        viewModelScope.launch {
            memoRepository.delete(memo)
        }
    }

    // Template management
    fun importTemplate(uri: Uri) {
        viewModelScope.launch {
            try {
                val context = getApplication<VoiceMemoApplication>()
                val inputStream = context.contentResolver.openInputStream(uri)
                val content = inputStream?.bufferedReader()?.readText() ?: return@launch
                inputStream.close()

                // Extract filename
                val name = uri.lastPathSegment
                    ?.substringAfterLast("/")
                    ?.removeSuffix(".md")
                    ?: "Geimporteerd sjabloon"

                val template = Template(
                    name = name,
                    content = content
                )
                templateRepository.insert(template)
            } catch (_: Exception) {
                // Handle error silently
            }
        }
    }

    fun saveTemplate(id: Long?, name: String, content: String, instructions: String) {
        viewModelScope.launch {
            if (id != null) {
                val existing = templateRepository.getTemplateById(id)
                if (existing != null) {
                    templateRepository.update(
                        existing.copy(name = name, content = content, instructions = instructions)
                    )
                }
            } else {
                templateRepository.insert(
                    Template(name = name, content = content, instructions = instructions)
                )
            }
        }
    }

    fun deleteTemplate(template: Template) {
        viewModelScope.launch {
            templateRepository.delete(template)
        }
    }

    override fun onCleared() {
        super.onCleared()
        speechService.destroy()
        geminiService.close()
    }
}
