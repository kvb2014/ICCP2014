package com.voicememo.app.ui.screens

import androidx.compose.animation.animateColorAsState
import androidx.compose.animation.core.*
import androidx.compose.foundation.background
import androidx.compose.foundation.layout.*
import androidx.compose.foundation.rememberScrollState
import androidx.compose.foundation.shape.CircleShape
import androidx.compose.foundation.verticalScroll
import androidx.compose.material.icons.Icons
import androidx.compose.material.icons.filled.Mic
import androidx.compose.material.icons.filled.Stop
import androidx.compose.material3.*
import androidx.compose.runtime.*
import androidx.compose.ui.Alignment
import androidx.compose.ui.Modifier
import androidx.compose.ui.draw.scale
import androidx.compose.ui.graphics.Color
import androidx.compose.ui.text.style.TextAlign
import androidx.compose.ui.unit.dp
import com.voicememo.app.service.SpeechRecognitionService.RecognitionState

@OptIn(ExperimentalMaterial3Api::class)
@Composable
fun RecordScreen(
    recognitionState: RecognitionState,
    currentTranscript: String,
    selectedTemplateName: String?,
    onStartRecording: () -> Unit,
    onStopRecording: () -> Unit,
    onSaveAndProcess: () -> Unit,
    onSelectTemplate: () -> Unit,
    modifier: Modifier = Modifier
) {
    val isListening = recognitionState is RecognitionState.Listening ||
            recognitionState is RecognitionState.Partial ||
            recognitionState is RecognitionState.Result

    val pulseAnimation = rememberInfiniteTransition(label = "pulse")
    val scale by pulseAnimation.animateFloat(
        initialValue = 1f,
        targetValue = 1.15f,
        animationSpec = infiniteRepeatable(
            animation = tween(800, easing = EaseInOutSine),
            repeatMode = RepeatMode.Reverse
        ),
        label = "scale"
    )

    val buttonColor by animateColorAsState(
        targetValue = if (isListening) MaterialTheme.colorScheme.error
        else MaterialTheme.colorScheme.primary,
        label = "buttonColor"
    )

    Scaffold(
        topBar = {
            TopAppBar(
                title = { Text("Spraakmemo") }
            )
        }
    ) { padding ->
        Column(
            modifier = modifier
                .fillMaxSize()
                .padding(padding)
                .padding(16.dp),
            horizontalAlignment = Alignment.CenterHorizontally
        ) {
            // Template selector
            OutlinedCard(
                onClick = onSelectTemplate,
                modifier = Modifier.fillMaxWidth()
            ) {
                Row(
                    modifier = Modifier
                        .fillMaxWidth()
                        .padding(16.dp),
                    horizontalArrangement = Arrangement.SpaceBetween,
                    verticalAlignment = Alignment.CenterVertically
                ) {
                    Column {
                        Text(
                            text = "Sjabloon",
                            style = MaterialTheme.typography.labelMedium,
                            color = MaterialTheme.colorScheme.onSurfaceVariant
                        )
                        Text(
                            text = selectedTemplateName ?: "Geen sjabloon geselecteerd",
                            style = MaterialTheme.typography.bodyLarge
                        )
                    }
                }
            }

            Spacer(modifier = Modifier.height(24.dp))

            // Transcript display
            Card(
                modifier = Modifier
                    .fillMaxWidth()
                    .weight(1f)
            ) {
                Column(
                    modifier = Modifier
                        .fillMaxSize()
                        .padding(16.dp)
                ) {
                    Text(
                        text = "Transcript",
                        style = MaterialTheme.typography.titleSmall,
                        color = MaterialTheme.colorScheme.onSurfaceVariant
                    )
                    Spacer(modifier = Modifier.height(8.dp))

                    val displayText = when (recognitionState) {
                        is RecognitionState.Partial -> recognitionState.text
                        is RecognitionState.Result -> recognitionState.text
                        else -> currentTranscript
                    }

                    if (displayText.isBlank()) {
                        Text(
                            text = if (isListening) "Luisteren..." else "Druk op de microfoon om te beginnen",
                            style = MaterialTheme.typography.bodyMedium,
                            color = MaterialTheme.colorScheme.onSurfaceVariant,
                            modifier = Modifier.fillMaxWidth(),
                            textAlign = TextAlign.Center
                        )
                    } else {
                        Text(
                            text = displayText,
                            style = MaterialTheme.typography.bodyMedium,
                            modifier = Modifier.verticalScroll(rememberScrollState())
                        )
                    }
                }
            }

            Spacer(modifier = Modifier.height(24.dp))

            // Status indicator
            if (recognitionState is RecognitionState.Error) {
                Text(
                    text = recognitionState.message,
                    color = MaterialTheme.colorScheme.error,
                    style = MaterialTheme.typography.bodySmall,
                    modifier = Modifier.padding(bottom = 8.dp)
                )
            }

            // Big microphone button
            Box(
                contentAlignment = Alignment.Center,
                modifier = Modifier.padding(bottom = 16.dp)
            ) {
                // Pulse ring when recording
                if (isListening) {
                    Box(
                        modifier = Modifier
                            .size(100.dp)
                            .scale(scale)
                            .background(
                                color = buttonColor.copy(alpha = 0.2f),
                                shape = CircleShape
                            )
                    )
                }

                FloatingActionButton(
                    onClick = {
                        if (isListening) onStopRecording() else onStartRecording()
                    },
                    modifier = Modifier.size(80.dp),
                    containerColor = buttonColor,
                    shape = CircleShape
                ) {
                    Icon(
                        imageVector = if (isListening) Icons.Default.Stop else Icons.Default.Mic,
                        contentDescription = if (isListening) "Stop opname" else "Start opname",
                        modifier = Modifier.size(40.dp),
                        tint = Color.White
                    )
                }
            }

            // Save button
            if (currentTranscript.isNotBlank() && !isListening) {
                Button(
                    onClick = onSaveAndProcess,
                    modifier = Modifier.fillMaxWidth()
                ) {
                    Text("Opslaan & Verwerken")
                }
            }
        }
    }
}
