package com.voicememo.app.data.db

import androidx.room.*
import com.voicememo.app.data.model.Memo
import kotlinx.coroutines.flow.Flow

@Dao
interface MemoDao {
    @Query("SELECT * FROM memos ORDER BY updatedAt DESC")
    fun getAllMemos(): Flow<List<Memo>>

    @Query("SELECT * FROM memos WHERE id = :id")
    suspend fun getMemoById(id: Long): Memo?

    @Insert
    suspend fun insert(memo: Memo): Long

    @Update
    suspend fun update(memo: Memo)

    @Delete
    suspend fun delete(memo: Memo)
}
