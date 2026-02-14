package com.voicememo.app.service

import android.os.Environment
import java.io.File

/**
 * Service responsible for reading/writing markdown files in the Obsidian vault.
 * Default vault path: /storage/emulated/0/Obsidian/
 */
class ObsidianExportService {

    companion object {
        val OBSIDIAN_BASE_PATH: String =
            Environment.getExternalStorageDirectory().absolutePath + "/Obsidian"
    }

    /**
     * Lists all vault directories under the Obsidian base path.
     */
    fun listVaults(): List<String> {
        val base = File(OBSIDIAN_BASE_PATH)
        if (!base.exists()) return emptyList()
        return base.listFiles()
            ?.filter { it.isDirectory && !it.name.startsWith(".") }
            ?.map { it.name }
            ?: emptyList()
    }

    /**
     * Saves markdown content to a file in the specified vault.
     * Returns the absolute path of the saved file.
     */
    fun saveNote(vaultName: String, fileName: String, content: String): String {
        val sanitizedName = sanitizeFileName(fileName)
        val vaultDir = File(OBSIDIAN_BASE_PATH, vaultName)
        if (!vaultDir.exists()) {
            vaultDir.mkdirs()
        }
        val file = File(vaultDir, "$sanitizedName.md")
        file.writeText(content)
        return file.absolutePath
    }

    /**
     * Reads an existing note from the vault.
     * Returns null if the file doesn't exist.
     */
    fun readNote(vaultName: String, fileName: String): String? {
        val sanitizedName = sanitizeFileName(fileName)
        val file = File(File(OBSIDIAN_BASE_PATH, vaultName), "$sanitizedName.md")
        return if (file.exists()) file.readText() else null
    }

    /**
     * Lists all markdown files in a vault.
     */
    fun listNotes(vaultName: String): List<String> {
        val vaultDir = File(OBSIDIAN_BASE_PATH, vaultName)
        if (!vaultDir.exists()) return emptyList()
        return vaultDir.walkTopDown()
            .filter { it.isFile && it.extension == "md" }
            .map { it.relativeTo(vaultDir).path.removeSuffix(".md") }
            .toList()
            .sorted()
    }

    /**
     * Checks if the Obsidian directory exists and is accessible.
     */
    fun isObsidianAvailable(): Boolean {
        return File(OBSIDIAN_BASE_PATH).exists()
    }

    private fun sanitizeFileName(name: String): String {
        return name.replace(Regex("[\\\\/:*?\"<>|]"), "_")
    }
}
