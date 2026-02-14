package com.voicememo.app.ui.screens

import androidx.compose.foundation.horizontalScroll
import androidx.compose.foundation.layout.*
import androidx.compose.foundation.rememberScrollState
import androidx.compose.foundation.verticalScroll
import androidx.compose.material.icons.Icons
import androidx.compose.material.icons.automirrored.filled.ArrowBack
import androidx.compose.material.icons.filled.Edit
import androidx.compose.material.icons.filled.Mic
import androidx.compose.material.icons.filled.Save
import androidx.compose.material3.*
import androidx.compose.runtime.*
import androidx.compose.ui.Modifier
import androidx.compose.ui.unit.dp
import com.voicememo.app.data.model.Memo

@OptIn(ExperimentalMaterial3Api::class)
@Composable
fun NoteDetailScreen(
    memo: Memo?,
    isProcessing: Boolean,
    onBack: () -> Unit,
    onContinueRecording: () -> Unit,
    onReprocess: () -> Unit,
    onEditMarkdown: (String) -> Unit,
    modifier: Modifier = Modifier
) {
    var isEditing by remember { mutableStateOf(false) }
    var editedMarkdown by remember(memo?.processedMarkdown) {
        mutableStateOf(memo?.processedMarkdown ?: "")
    }
    var selectedTab by remember { mutableIntStateOf(0) }

    Scaffold(
        topBar = {
            TopAppBar(
                title = { Text(memo?.title ?: "Laden...") },
                navigationIcon = {
                    IconButton(onClick = onBack) {
                        Icon(Icons.AutoMirrored.Filled.ArrowBack, "Terug")
                    }
                },
                actions = {
                    IconButton(onClick = onContinueRecording) {
                        Icon(Icons.Default.Mic, "Verder opnemen")
                    }
                    if (selectedTab == 1) {
                        if (isEditing) {
                            IconButton(onClick = {
                                isEditing = false
                                onEditMarkdown(editedMarkdown)
                            }) {
                                Icon(Icons.Default.Save, "Opslaan")
                            }
                        } else {
                            IconButton(onClick = { isEditing = true }) {
                                Icon(Icons.Default.Edit, "Bewerken")
                            }
                        }
                    }
                }
            )
        }
    ) { padding ->
        Column(
            modifier = modifier
                .fillMaxSize()
                .padding(padding)
        ) {
            TabRow(selectedTabIndex = selectedTab) {
                Tab(
                    selected = selectedTab == 0,
                    onClick = { selectedTab = 0 },
                    text = { Text("Transcript") }
                )
                Tab(
                    selected = selectedTab == 1,
                    onClick = { selectedTab = 1 },
                    text = { Text("Aantekening") }
                )
            }

            when (selectedTab) {
                0 -> {
                    // Transcript tab
                    Column(
                        modifier = Modifier
                            .fillMaxSize()
                            .padding(16.dp)
                            .verticalScroll(rememberScrollState())
                    ) {
                        Text(
                            text = memo?.transcript ?: "",
                            style = MaterialTheme.typography.bodyMedium
                        )
                    }
                }

                1 -> {
                    // Markdown note tab
                    if (isProcessing) {
                        Box(
                            modifier = Modifier.fillMaxSize(),
                            contentAlignment = androidx.compose.ui.Alignment.Center
                        ) {
                            Column(
                                horizontalAlignment = androidx.compose.ui.Alignment.CenterHorizontally
                            ) {
                                CircularProgressIndicator()
                                Spacer(modifier = Modifier.height(16.dp))
                                Text("Verwerken met Gemini Nano...")
                            }
                        }
                    } else if (memo?.processedMarkdown.isNullOrBlank()) {
                        Box(
                            modifier = Modifier.fillMaxSize(),
                            contentAlignment = androidx.compose.ui.Alignment.Center
                        ) {
                            Column(
                                horizontalAlignment = androidx.compose.ui.Alignment.CenterHorizontally
                            ) {
                                Text(
                                    text = "Nog niet verwerkt",
                                    style = MaterialTheme.typography.bodyLarge
                                )
                                Spacer(modifier = Modifier.height(8.dp))
                                Button(onClick = onReprocess) {
                                    Text("Nu verwerken")
                                }
                            }
                        }
                    } else if (isEditing) {
                        OutlinedTextField(
                            value = editedMarkdown,
                            onValueChange = { editedMarkdown = it },
                            modifier = Modifier
                                .fillMaxSize()
                                .padding(8.dp),
                            textStyle = MaterialTheme.typography.bodySmall
                        )
                    } else {
                        Column(
                            modifier = Modifier
                                .fillMaxSize()
                                .padding(16.dp)
                                .verticalScroll(rememberScrollState())
                        ) {
                            Text(
                                text = memo?.processedMarkdown ?: "",
                                style = MaterialTheme.typography.bodyMedium
                            )
                        }
                    }
                }
            }
        }
    }
}
