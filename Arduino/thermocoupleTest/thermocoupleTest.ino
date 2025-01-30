#include "max6675.h"

int thermoDO = 12;
int thermoCS = 37;
int thermoCLK = 13;
float T = 0;

MAX6675 thermocouple(thermoCLK, thermoCS, thermoDO);

void setup() {
  Serial.begin(9600);

  Serial.println("MAX6675 test");
  // wait for MAX chip to stabilize
  delay(500);
}

void loop() {
  // basic readout test, just print the current temp
  
  T = thermocouple.readFahrenheit();
  Serial.println(thermocouple.readFahrenheit());



/*  T_air = thermocouple.readFahrenheit();
  offset = T_boiling - T_air
  T_calibrated = T_uncalibrated
*/
 
  // For the MAX6675 to update, you must delay AT LEAST 250ms between reads!
  delay(250);
}