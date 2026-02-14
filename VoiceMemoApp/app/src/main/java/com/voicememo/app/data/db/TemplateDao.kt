package com.voicememo.app.data.db

import androidx.room.*
import com.voicememo.app.data.model.Template
import kotlinx.coroutines.flow.Flow

@Dao
interface TemplateDao {
    @Query("SELECT * FROM templates ORDER BY name ASC")
    fun getAllTemplates(): Flow<List<Template>>

    @Query("SELECT * FROM templates WHERE id = :id")
    suspend fun getTemplateById(id: Long): Template?

    @Insert
    suspend fun insert(template: Template): Long

    @Update
    suspend fun update(template: Template)

    @Delete
    suspend fun delete(template: Template)
}
