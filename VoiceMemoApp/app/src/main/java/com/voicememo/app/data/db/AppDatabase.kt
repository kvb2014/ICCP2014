package com.voicememo.app.data.db

import android.content.Context
import androidx.room.Database
import androidx.room.Room
import androidx.room.RoomDatabase
import com.voicememo.app.data.model.Memo
import com.voicememo.app.data.model.Template

@Database(entities = [Memo::class, Template::class], version = 1, exportSchema = false)
abstract class AppDatabase : RoomDatabase() {
    abstract fun memoDao(): MemoDao
    abstract fun templateDao(): TemplateDao

    companion object {
        @Volatile
        private var INSTANCE: AppDatabase? = null

        fun getDatabase(context: Context): AppDatabase {
            return INSTANCE ?: synchronized(this) {
                val instance = Room.databaseBuilder(
                    context.applicationContext,
                    AppDatabase::class.java,
                    "voicememo_database"
                ).build()
                INSTANCE = instance
                instance
            }
        }
    }
}
