package com.voicememo.app.data.model

import androidx.room.Entity
import androidx.room.PrimaryKey

@Entity(tableName = "memos")
data class Memo(
    @PrimaryKey(autoGenerate = true)
    val id: Long = 0,
    val title: String,
    val transcript: String,
    val processedMarkdown: String? = null,
    val templateId: Long? = null,
    val obsidianPath: String? = null,
    val createdAt: Long = System.currentTimeMillis(),
    val updatedAt: Long = System.currentTimeMillis()
)
