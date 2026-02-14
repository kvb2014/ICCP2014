package com.voicememo.app.ui.screens

import android.net.Uri
import androidx.activity.compose.rememberLauncherForActivityResult
import androidx.activity.result.contract.ActivityResultContracts
import androidx.compose.foundation.layout.*
import androidx.compose.foundation.lazy.LazyColumn
import androidx.compose.foundation.lazy.items
import androidx.compose.material.icons.Icons
import androidx.compose.material.icons.filled.Add
import androidx.compose.material.icons.filled.Delete
import androidx.compose.material.icons.filled.Edit
import androidx.compose.material.icons.filled.UploadFile
import androidx.compose.material3.*
import androidx.compose.runtime.*
import androidx.compose.ui.Alignment
import androidx.compose.ui.Modifier
import androidx.compose.ui.text.style.TextOverflow
import androidx.compose.ui.unit.dp
import com.voicememo.app.data.model.Template

@OptIn(ExperimentalMaterial3Api::class)
@Composable
fun TemplateListScreen(
    templates: List<Template>,
    onImportFile: (Uri) -> Unit,
    onCreateTemplate: () -> Unit,
    onEditTemplate: (Template) -> Unit,
    onDeleteTemplate: (Template) -> Unit,
    modifier: Modifier = Modifier
) {
    val filePickerLauncher = rememberLauncherForActivityResult(
        contract = ActivityResultContracts.OpenDocument()
    ) { uri: Uri? ->
        uri?.let { onImportFile(it) }
    }

    Scaffold(
        topBar = {
            TopAppBar(
                title = { Text("Sjablonen") }
            )
        },
        floatingActionButton = {
            Column(
                horizontalAlignment = Alignment.End,
                verticalArrangement = Arrangement.spacedBy(12.dp)
            ) {
                SmallFloatingActionButton(
                    onClick = { filePickerLauncher.launch(arrayOf("text/*", "*/*")) }
                ) {
                    Icon(Icons.Default.UploadFile, "Importeer .md bestand")
                }
                FloatingActionButton(onClick = onCreateTemplate) {
                    Icon(Icons.Default.Add, "Nieuw sjabloon")
                }
            }
        }
    ) { padding ->
        if (templates.isEmpty()) {
            Box(
                modifier = modifier
                    .fillMaxSize()
                    .padding(padding),
                contentAlignment = Alignment.Center
            ) {
                Column(horizontalAlignment = Alignment.CenterHorizontally) {
                    Text(
                        text = "Geen sjablonen",
                        style = MaterialTheme.typography.headlineSmall,
                        color = MaterialTheme.colorScheme.onSurfaceVariant
                    )
                    Spacer(modifier = Modifier.height(8.dp))
                    Text(
                        text = "Upload een .md bestand of maak een nieuw sjabloon aan",
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
                items(templates, key = { it.id }) { template ->
                    TemplateCard(
                        template = template,
                        onEdit = { onEditTemplate(template) },
                        onDelete = { onDeleteTemplate(template) }
                    )
                }
            }
        }
    }
}

@Composable
private fun TemplateCard(
    template: Template,
    onEdit: () -> Unit,
    onDelete: () -> Unit
) {
    var showDeleteDialog by remember { mutableStateOf(false) }

    Card(
        modifier = Modifier.fillMaxWidth()
    ) {
        Column(modifier = Modifier.padding(16.dp)) {
            Row(
                modifier = Modifier.fillMaxWidth(),
                horizontalArrangement = Arrangement.SpaceBetween,
                verticalAlignment = Alignment.CenterVertically
            ) {
                Text(
                    text = template.name,
                    style = MaterialTheme.typography.titleMedium,
                    modifier = Modifier.weight(1f)
                )
                Row {
                    IconButton(onClick = onEdit) {
                        Icon(Icons.Default.Edit, "Bewerk")
                    }
                    IconButton(onClick = { showDeleteDialog = true }) {
                        Icon(Icons.Default.Delete, "Verwijder")
                    }
                }
            }
            if (template.instructions.isNotBlank()) {
                Spacer(modifier = Modifier.height(4.dp))
                Text(
                    text = template.instructions,
                    style = MaterialTheme.typography.bodySmall,
                    color = MaterialTheme.colorScheme.onSurfaceVariant,
                    maxLines = 2,
                    overflow = TextOverflow.Ellipsis
                )
            }
            Spacer(modifier = Modifier.height(4.dp))
            Text(
                text = template.content.take(150) + if (template.content.length > 150) "..." else "",
                style = MaterialTheme.typography.bodySmall,
                color = MaterialTheme.colorScheme.onSurfaceVariant,
                maxLines = 3,
                overflow = TextOverflow.Ellipsis
            )
        }
    }

    if (showDeleteDialog) {
        AlertDialog(
            onDismissRequest = { showDeleteDialog = false },
            title = { Text("Sjabloon verwijderen?") },
            text = { Text("Weet je zeker dat je '${template.name}' wilt verwijderen?") },
            confirmButton = {
                TextButton(onClick = {
                    showDeleteDialog = false
                    onDelete()
                }) {
                    Text("Verwijderen")
                }
            },
            dismissButton = {
                TextButton(onClick = { showDeleteDialog = false }) {
                    Text("Annuleren")
                }
            }
        )
    }
}
