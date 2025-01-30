int analogPin = A3;
int val = 0;
int refVoltage = 5; // 5VDC
int numBits = 10;
int res = pow(2,numBits) - 1;
float pressure;

void setup() {
  // put your setup code here, to run once:
  Serial.begin(9600);
}

void loop() {
  // put your main code here, to run repeatedly:
  
  // get digital value corresponding
  val = analogRead(analogPin);
  //Serial.println(val);
  //Serial.println(res);

  pressure = val * 1000.0 / 1024;
  Serial.println(pressure);



}
