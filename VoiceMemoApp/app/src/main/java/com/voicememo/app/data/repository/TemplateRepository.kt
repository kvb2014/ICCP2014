package com.voicememo.app.data.repository

import com.voicememo.app.data.db.TemplateDao
import com.voicememo.app.data.model.Template
import kotlinx.coroutines.flow.Flow

class TemplateRepository(private val templateDao: TemplateDao) {
    val allTemplates: Flow<List<Template>> = templateDao.getAllTemplates()

    suspend fun getTemplateById(id: Long): Template? = templateDao.getTemplateById(id)

    suspend fun insert(template: Template): Long = templateDao.insert(template)

    suspend fun update(template: Template) = templateDao.update(template)

    suspend fun delete(template: Template) = templateDao.delete(template)
}
