package com.voicememo.app.data.repository

import com.voicememo.app.data.db.MemoDao
import com.voicememo.app.data.model.Memo
import kotlinx.coroutines.flow.Flow

class MemoRepository(private val memoDao: MemoDao) {
    val allMemos: Flow<List<Memo>> = memoDao.getAllMemos()

    suspend fun getMemoById(id: Long): Memo? = memoDao.getMemoById(id)

    suspend fun insert(memo: Memo): Long = memoDao.insert(memo)

    suspend fun update(memo: Memo) = memoDao.update(memo)

    suspend fun delete(memo: Memo) = memoDao.delete(memo)
}
