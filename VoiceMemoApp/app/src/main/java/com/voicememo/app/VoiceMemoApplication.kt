package com.voicememo.app

import android.app.Application
import com.voicememo.app.data.db.AppDatabase
import com.voicememo.app.data.repository.MemoRepository
import com.voicememo.app.data.repository.TemplateRepository

class VoiceMemoApplication : Application() {
    val database by lazy { AppDatabase.getDatabase(this) }
    val memoRepository by lazy { MemoRepository(database.memoDao()) }
    val templateRepository by lazy { TemplateRepository(database.templateDao()) }
}
