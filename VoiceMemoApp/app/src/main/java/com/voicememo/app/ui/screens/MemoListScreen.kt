package com.voicememo.app.ui.screens

import androidx.compose.foundation.layout.*
import androidx.compose.foundation.lazy.LazyColumn
import androidx.compose.foundation.lazy.items
import androidx.compose.material.icons.Icons
import androidx.compose.material.icons.filled.Delete
import androidx.compose.material.icons.filled.Mic
import androidx.compose.material.icons.filled.NoteAlt
import androidx.compose.material3.*
import androidx.compose.runtime.*
import androidx.compose.ui.Alignment
import androidx.compose.ui.Modifier
import androidx.compose.ui.text.style.TextOverflow
import androidx.compose.ui.unit.dp
import com.voicememo.app.data.model.Memo
import java.text.SimpleDateFormat
import java.util.*

@OptIn(ExperimentalMaterial3Api::class)
@Composable
fun MemoListScreen(
    memos: List<Memo>,
    onMemoClick: (Memo) -> Unit,
    onContinueRecording: (Memo) -> Unit,
    onDeleteMemo: (Memo) -> Unit,
    modifier: Modifier = Modifier
) {
    Scaffold(
        topBar = {
            TopAppBar(
                title = { Text("Memo's") }
            )
        }
    ) { padding ->
        if (memos.isEmpty()) {
            Box(
                modifier = modifier
                    .fillMaxSize()
                    .padding(padding),
                contentAlignment = Alignment.Center
            ) {
                Column(horizontalAlignment = Alignment.CenterHorizontally) {
                    Text(
                        text = "Geen memo's",
                        style = MaterialTheme.typography.headlineSmall,
                        color = MaterialTheme.colorScheme.onSurfaceVariant
                    )
                    Spacer(modifier = Modifier.height(8.dp))
                    Text(
                        text = "Neem een spraakmemo op om te beginnen",
                        style = MaterialTheme.typography.bodyMedium,
                        color = MaterialTheme.colorScheme.onSurfaceVariant
                    )
                }
            }
        } else {
            LazyColumn(
                modifier = modifier
                    .fillMaxSize()
                    .padding(padding),
                contentPadding = PaddingValues(16.dp),
                verticalArrangement = Arrangement.spacedBy(8.dp)
            ) {
                items(memos, key = { it.id }) { memo ->
                    MemoCard(
                        memo = memo,
                        onClick = { onMemoClick(memo) },
                        onContinue = { onContinueRecording(memo) },
                        onDelete = { onDeleteMemo(memo) }
                    )
                }
            }
        }
    }
}

@Composable
private fun MemoCard(
    memo: Memo,
    onClick: () -> Unit,
    onContinue: () -> Unit,
    onDelete: () -> Unit
) {
    var showDeleteDialog by remember { mutableStateOf(false) }
    val dateFormat = remember { SimpleDateFormat("dd-MM-yyyy HH:mm", Locale("nl")) }

    Card(
        onClick = onClick,
        modifier = Modifier.fillMaxWidth()
    ) {
        Column(modifier = Modifier.padding(16.dp)) {
            Row(
                modifier = Modifier.fillMaxWidth(),
                horizontalArrangement = Arrangement.SpaceBetween,
                verticalAlignment = Alignment.CenterVertically
            ) {
                Column(modifier = Modifier.weight(1f)) {
                    Text(
                        text = memo.title,
                        style = MaterialTheme.typography.titleMedium
                    )
                    Text(
                        text = dateFormat.format(Date(memo.updatedAt)),
                        style = MaterialTheme.typography.bodySmall,
                        color = MaterialTheme.colorScheme.onSurfaceVariant
                    )
                }
                Row {
                    IconButton(onClick = onContinue) {
                        Icon(
                            Icons.Default.Mic,
                            contentDescription = "Verder opnemen",
                            tint = MaterialTheme.colorScheme.primary
                        )
                    }
                    IconButton(onClick = { showDeleteDialog = true }) {
                        Icon(Icons.Default.Delete, "Verwijder")
                    }
                }
            }
            Spacer(modifier = Modifier.height(8.dp))
            Text(
                text = memo.transcript.take(120) + if (memo.transcript.length > 120) "..." else "",
                style = MaterialTheme.typography.bodySmall,
                color = MaterialTheme.colorScheme.onSurfaceVariant,
                maxLines = 2,
                overflow = TextOverflow.Ellipsis
            )
            if (memo.processedMarkdown != null) {
                Spacer(modifier = Modifier.height(4.dp))
                Row(verticalAlignment = Alignment.CenterVertically) {
                    Icon(
                        Icons.Default.NoteAlt,
                        contentDescription = null,
                        modifier = Modifier.size(14.dp),
                        tint = MaterialTheme.colorScheme.primary
                    )
                    Spacer(modifier = Modifier.width(4.dp))
                    Text(
                        text = "Verwerkt naar Obsidian",
                        style = MaterialTheme.typography.labelSmall,
                        color = MaterialTheme.colorScheme.primary
                    )
                }
            }
        }
    }

    if (showDeleteDialog) {
        AlertDialog(
            onDismissRequest = { showDeleteDialog = false },
            title = { Text("Memo verwijderen?") },
            text = { Text("Weet je zeker dat je '${memo.title}' wilt verwijderen?") },
            confirmButton = {
                TextButton(onClick = {
                    showDeleteDialog = false
                    onDelete()
                }) { Text("Verwijderen") }
            },
            dismissButton = {
                TextButton(onClick = { showDeleteDialog = false }) { Text("Annuleren") }
            }
        )
    }
}
