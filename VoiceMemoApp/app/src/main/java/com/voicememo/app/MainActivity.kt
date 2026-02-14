package com.voicememo.app

import android.Manifest
import android.content.Intent
import android.content.pm.PackageManager
import android.net.Uri
import android.os.Build
import android.os.Bundle
import android.os.Environment
import android.provider.Settings
import androidx.activity.ComponentActivity
import androidx.activity.compose.setContent
import androidx.activity.result.contract.ActivityResultContracts
import androidx.compose.foundation.layout.*
import androidx.compose.material.icons.Icons
import androidx.compose.material.icons.filled.Description
import androidx.compose.material.icons.filled.List
import androidx.compose.material.icons.filled.Mic
import androidx.compose.material3.*
import androidx.compose.runtime.*
import androidx.compose.ui.Modifier
import androidx.core.content.ContextCompat
import androidx.lifecycle.viewmodel.compose.viewModel
import androidx.navigation.NavType
import androidx.navigation.compose.*
import androidx.navigation.navArgument
import com.voicememo.app.ui.MainViewModel
import com.voicememo.app.ui.screens.*
import com.voicememo.app.ui.theme.VoiceMemoTheme

class MainActivity : ComponentActivity() {

    private val requestPermissions = registerForActivityResult(
        ActivityResultContracts.RequestMultiplePermissions()
    ) { /* permissions handled */ }

    override fun onCreate(savedInstanceState: Bundle?) {
        super.onCreate(savedInstanceState)
        requestNecessaryPermissions()

        setContent {
            VoiceMemoTheme {
                VoiceMemoApp()
            }
        }
    }

    private fun requestNecessaryPermissions() {
        val permissions = mutableListOf(Manifest.permission.RECORD_AUDIO)

        if (Build.VERSION.SDK_INT < Build.VERSION_CODES.R) {
            permissions.add(Manifest.permission.WRITE_EXTERNAL_STORAGE)
            permissions.add(Manifest.permission.READ_EXTERNAL_STORAGE)
        }

        val needsPermission = permissions.any {
            ContextCompat.checkSelfPermission(this, it) != PackageManager.PERMISSION_GRANTED
        }

        if (needsPermission) {
            requestPermissions.launch(permissions.toTypedArray())
        }

        // For Android 11+, request MANAGE_EXTERNAL_STORAGE for Obsidian vault access
        if (Build.VERSION.SDK_INT >= Build.VERSION_CODES.R &&
            !Environment.isExternalStorageManager()
        ) {
            val intent = Intent(Settings.ACTION_MANAGE_APP_ALL_FILES_ACCESS_PERMISSION).apply {
                data = Uri.parse("package:$packageName")
            }
            startActivity(intent)
        }
    }
}

@Composable
fun VoiceMemoApp(viewModel: MainViewModel = viewModel()) {
    val navController = rememberNavController()
    val navBackStackEntry by navController.currentBackStackEntryAsState()
    val currentRoute = navBackStackEntry?.destination?.route

    // Dialogs
    var showTemplateSelector by remember { mutableStateOf(false) }
    var showVaultSelector by remember { mutableStateOf(false) }

    // State
    val recognitionState by viewModel.speechService.state.collectAsState()
    val currentTranscript by viewModel.currentTranscript.collectAsState()
    val templates by viewModel.templates.collectAsState()
    val memos by viewModel.memos.collectAsState()
    val selectedTemplate by viewModel.selectedTemplate.collectAsState()
    val selectedVault by viewModel.selectedVault.collectAsState()
    val isProcessing by viewModel.isProcessing.collectAsState()
    val detailMemo by viewModel.detailMemo.collectAsState()

    // Editing state for template edit screen
    var editingTemplate by remember { mutableStateOf<com.voicememo.app.data.model.Template?>(null) }
    var showTemplateEdit by remember { mutableStateOf(false) }

    Scaffold(
        bottomBar = {
            if (currentRoute in listOf("record", "memos", "templates")) {
                NavigationBar {
                    NavigationBarItem(
                        icon = { Icon(Icons.Default.Mic, contentDescription = "Opnemen") },
                        label = { Text("Opnemen") },
                        selected = currentRoute == "record",
                        onClick = {
                            navController.navigate("record") {
                                popUpTo("record") { inclusive = true }
                            }
                        }
                    )
                    NavigationBarItem(
                        icon = { Icon(Icons.Default.List, contentDescription = "Memo's") },
                        label = { Text("Memo's") },
                        selected = currentRoute == "memos",
                        onClick = {
                            navController.navigate("memos") {
                                popUpTo("record")
                            }
                        }
                    )
                    NavigationBarItem(
                        icon = {
                            Icon(
                                Icons.Default.Description,
                                contentDescription = "Sjablonen"
                            )
                        },
                        label = { Text("Sjablonen") },
                        selected = currentRoute == "templates",
                        onClick = {
                            navController.navigate("templates") {
                                popUpTo("record")
                            }
                        }
                    )
                }
            }
        }
    ) { innerPadding ->
        NavHost(
            navController = navController,
            startDestination = "record",
            modifier = Modifier.padding(innerPadding)
        ) {
            composable("record") {
                RecordScreen(
                    recognitionState = recognitionState,
                    currentTranscript = currentTranscript,
                    selectedTemplateName = selectedTemplate?.name,
                    onStartRecording = {
                        if (selectedVault == null) {
                            showVaultSelector = true
                        } else {
                            viewModel.startRecording()
                        }
                    },
                    onStopRecording = { viewModel.stopRecording() },
                    onSaveAndProcess = {
                        viewModel.saveAndProcess()
                        navController.navigate("memos")
                    },
                    onSelectTemplate = { showTemplateSelector = true }
                )
            }

            composable("memos") {
                MemoListScreen(
                    memos = memos,
                    onMemoClick = { memo ->
                        viewModel.loadMemoDetail(memo.id)
                        navController.navigate("detail/${memo.id}")
                    },
                    onContinueRecording = { memo ->
                        viewModel.continueRecording(memo)
                        navController.navigate("record")
                    },
                    onDeleteMemo = { viewModel.deleteMemo(it) }
                )
            }

            composable("templates") {
                TemplateListScreen(
                    templates = templates,
                    onImportFile = { uri -> viewModel.importTemplate(uri) },
                    onCreateTemplate = {
                        editingTemplate = null
                        showTemplateEdit = true
                        navController.navigate("template_edit")
                    },
                    onEditTemplate = { template ->
                        editingTemplate = template
                        showTemplateEdit = true
                        navController.navigate("template_edit")
                    },
                    onDeleteTemplate = { viewModel.deleteTemplate(it) }
                )
            }

            composable("template_edit") {
                TemplateEditScreen(
                    template = editingTemplate,
                    onSave = { name, content, instructions ->
                        viewModel.saveTemplate(
                            editingTemplate?.id,
                            name,
                            content,
                            instructions
                        )
                        navController.popBackStack()
                    },
                    onBack = { navController.popBackStack() }
                )
            }

            composable(
                "detail/{memoId}",
                arguments = listOf(navArgument("memoId") { type = NavType.LongType })
            ) { backStackEntry ->
                val memoId = backStackEntry.arguments?.getLong("memoId") ?: return@composable

                LaunchedEffect(memoId) {
                    viewModel.loadMemoDetail(memoId)
                }

                NoteDetailScreen(
                    memo = detailMemo,
                    isProcessing = isProcessing,
                    onBack = { navController.popBackStack() },
                    onContinueRecording = {
                        detailMemo?.let { memo ->
                            viewModel.continueRecording(memo)
                            navController.navigate("record") {
                                popUpTo("record") { inclusive = true }
                            }
                        }
                    },
                    onReprocess = { viewModel.reprocessMemo(memoId) },
                    onEditMarkdown = { markdown ->
                        viewModel.updateMemoMarkdown(memoId, markdown)
                    }
                )
            }
        }
    }

    // Template selector dialog
    if (showTemplateSelector) {
        TemplateSelectorDialog(
            templates = templates,
            onSelect = { template ->
                viewModel.selectTemplate(template)
                showTemplateSelector = false
            },
            onDismiss = { showTemplateSelector = false }
        )
    }

    // Vault selector dialog
    if (showVaultSelector) {
        VaultSelectorDialog(
            vaults = viewModel.obsidianService.listVaults(),
            onSelect = { vault ->
                viewModel.selectVault(vault)
                showVaultSelector = false
                viewModel.startRecording()
            },
            onDismiss = { showVaultSelector = false }
        )
    }
}
