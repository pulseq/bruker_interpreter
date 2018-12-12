# Pulseq to Bruker PPG converter

The code in this repository converts Pulseq .seq files (version 1.0.0) to Bruker .ppg files. 
Note that this may not work for you without modifications.

Due to the lack of manpower and time, we decided release this code to the community in the hope that others will continue the development.
Please consider issuing a pull request if you modify the code in a way which may be useful to others.

Ideally this converter should be replaced by a Paravision Method (PVM) which reads the pulseq file, converts it to PPG and WAVE files and sets all necessary variables.

## Disclaimer

We do not take any responsibility and we are not liable for any damage caused through the use of this software.

## Known problems

* No support for gradient shapes yet
* Long sequences (for instance 3D) can create PPG files too large for the scanner

## License

MIT License ([LICENSE](LICENSE) or http://opensource.org/licenses/MIT)
