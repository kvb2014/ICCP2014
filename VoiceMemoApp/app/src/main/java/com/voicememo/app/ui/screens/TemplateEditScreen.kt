package com.voicememo.app.ui.screens

import androidx.compose.foundation.layout.*
import androidx.compose.foundation.rememberScrollState
import androidx.compose.foundation.verticalScroll
import androidx.compose.material.icons.Icons
import androidx.compose.material.icons.automirrored.filled.ArrowBack
import androidx.compose.material.icons.filled.Save
import androidx.compose.material3.*
import androidx.compose.runtime.*
import androidx.compose.ui.Modifier
import androidx.compose.ui.unit.dp
import com.voicememo.app.data.model.Template

@OptIn(ExperimentalMaterial3Api::class)
@Composable
fun TemplateEditScreen(
    template: Template?,
    onSave: (name: String, content: String, instructions: String) -> Unit,
    onBack: () -> Unit,
    modifier: Modifier = Modifier
) {
    var name by remember { mutableStateOf(template?.name ?: "") }
    var content by remember { mutableStateOf(template?.content ?: "") }
    var instructions by remember { mutableStateOf(template?.instructions ?: "") }

    Scaffold(
        topBar = {
            TopAppBar(
                title = { Text(if (template == null) "Nieuw sjabloon" else "Sjabloon bewerken") },
                navigationIcon = {
                    IconButton(onClick = onBack) {
                        Icon(Icons.AutoMirrored.Filled.ArrowBack, "Terug")
                    }
                },
                actions = {
                    IconButton(
                        onClick = { onSave(name, content, instructions) },
                        enabled = name.isNotBlank() && content.isNotBlank()
                    ) {
                        Icon(Icons.Default.Save, "Opslaan")
                    }
                }
            )
        }
    ) { padding ->
        Column(
            modifier = modifier
                .fillMaxSize()
                .padding(padding)
                .padding(16.dp)
                .verticalScroll(rememberScrollState()),
            verticalArrangement = Arrangement.spacedBy(16.dp)
        ) {
            OutlinedTextField(
                value = name,
                onValueChange = { name = it },
                label = { Text("Naam") },
                placeholder = { Text("Bijv. Podcast notities") },
                singleLine = true,
                modifier = Modifier.fillMaxWidth()
            )

            OutlinedTextField(
                value = instructions,
                onValueChange = { instructions = it },
                label = { Text("Instructies voor de AI") },
                placeholder = {
                    Text(
                        "Bijv. Vat de belangrijkste punten samen. " +
                                "Gebruik tussenkopjes voor de hoofdonderwerpen."
                    )
                },
                minLines = 3,
                maxLines = 6,
                modifier = Modifier.fillMaxWidth()
            )

            OutlinedTextField(
                value = content,
                onValueChange = { content = it },
                label = { Text("Markdown sjabloon") },
                placeholder = {
                    Text(
                        "---\ntags: [notitie]\ndatum: {{datum}}\n---\n# {{titel}}\n\n## Samenvatting\n..."
                    )
                },
                minLines = 10,
                maxLines = 20,
                modifier = Modifier.fillMaxWidth()
            )

            Text(
                text = "Tip: Gebruik {{titel}}, {{datum}}, {{transcript}} als placeholders. " +
                        "De AI vult de rest intelligent in op basis van het transcript.",
                style = MaterialTheme.typography.bodySmall,
                color = MaterialTheme.colorScheme.onSurfaceVariant
            )
        }
    }
}
