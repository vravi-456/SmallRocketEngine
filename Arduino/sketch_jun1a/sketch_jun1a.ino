int relayPin = 11;
boolean run = true;

void setup() {
  // put your setup code here, to run once:
  pinMode(relayPin, OUTPUT);
}

void loop() {
  // put your main code here, to run repeatedly:
  digitalWrite(relayPin, LOW);
  delay(5000);
  Serial.print("Here");
}
