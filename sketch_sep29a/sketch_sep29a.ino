// Configuración para Arduino
const int analogPin = A0;
const int sampleRate = 1000; // Frecuencia de muestreo (Hz)
const int numSamples = 5000; // Número de muestras a capturar

void setup() {
  Serial.begin(115200); // Inicializar comunicación serie
}

void loop() {
  for (int i = 0; i < numSamples; i++) {
    int analogValue = analogRead(analogPin); // Leer pin A0
    Serial.println(analogValue); // Enviar dato a MATLAB
    delayMicroseconds(1000000 / sampleRate); // Control de la tasa de muestreo
  }
  delay(1000); // Pausa entre lecturas
}

