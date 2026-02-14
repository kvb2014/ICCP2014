package com.voicememo.app.ui.screens

import androidx.compose.foundation.clickable
import androidx.compose.foundation.layout.*
import androidx.compose.foundation.lazy.LazyColumn
import androidx.compose.foundation.lazy.items
import androidx.compose.material3.*
import androidx.compose.runtime.Composable
import androidx.compose.ui.Modifier

@Composable
fun VaultSelectorDialog(
    vaults: List<String>,
    onSelect: (String) -> Unit,
    onDismiss: () -> Unit
) {
    AlertDialog(
        onDismissRequest = onDismiss,
        title = { Text("Kies een Obsidian vault") },
        text = {
            if (vaults.isEmpty()) {
                Text("Geen Obsidian vaults gevonden in /storage/emulated/0/Obsidian/. Zorg dat Obsidian is geinstalleerd en een vault is aangemaakt.")
            } else {
                LazyColumn {
                    items(vaults) { vault ->
                        ListItem(
                            headlineContent = { Text(vault) },
                            modifier = Modifier.clickable { onSelect(vault) }
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
