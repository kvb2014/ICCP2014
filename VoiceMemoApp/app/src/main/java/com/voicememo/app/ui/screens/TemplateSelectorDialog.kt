package com.voicememo.app.ui.screens

import androidx.compose.foundation.clickable
import androidx.compose.foundation.layout.*
import androidx.compose.foundation.lazy.LazyColumn
import androidx.compose.foundation.lazy.items
import androidx.compose.material3.*
import androidx.compose.runtime.Composable
import androidx.compose.ui.Modifier
import androidx.compose.ui.unit.dp
import com.voicememo.app.data.model.Template

@Composable
fun TemplateSelectorDialog(
    templates: List<Template>,
    onSelect: (Template) -> Unit,
    onDismiss: () -> Unit
) {
    AlertDialog(
        onDismissRequest = onDismiss,
        title = { Text("Kies een sjabloon") },
        text = {
            if (templates.isEmpty()) {
                Text("Geen sjablonen beschikbaar. Maak er eerst een aan via het Sjablonen-tabblad.")
            } else {
                LazyColumn {
                    items(templates) { template ->
                        ListItem(
                            headlineContent = { Text(template.name) },
                            supportingContent = {
                                if (template.instructions.isNotBlank()) {
                                    Text(
                                        template.instructions.take(80),
                                        maxLines = 1
                                    )
                                }
                            },
                            modifier = Modifier.clickable { onSelect(template) }
                        )
                        HorizontalDivider()
                    }
                }
            }
        },
        confirmButton = {
            TextButton(onClick = onDismiss) { Text("Annuleren") }
        }
    )
}
