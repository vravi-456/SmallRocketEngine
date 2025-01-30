int relayPin = 13;

void setup() {
  // put your setup code here, to run once:
  pinMode(relayPin, OUTPUT);
}

void loop() {
  // put your main code here, to run repeatedly:

  // close for 3 seconds
  digitalWrite(relayPin, HIGH);
  delay(3000);

  // open for 1000 seconds
  digitalWrite(relayPin, LOW);
  delay(1000000);

  // close for two seconds
  digitalWrite(relayPin, HIGH);
  delay(2000);
}
