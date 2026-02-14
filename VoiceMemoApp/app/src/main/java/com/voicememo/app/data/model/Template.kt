package com.voicememo.app.data.model

import androidx.room.Entity
import androidx.room.PrimaryKey

@Entity(tableName = "templates")
data class Template(
    @PrimaryKey(autoGenerate = true)
    val id: Long = 0,
    val name: String,
    val content: String,
    val instructions: String = "",
    val createdAt: Long = System.currentTimeMillis()
)
