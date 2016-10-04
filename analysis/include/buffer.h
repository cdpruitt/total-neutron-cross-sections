// Raw data is stored as hexadecimal words (16 bits long) in the .evt files

// to read from a filestream, we need a buffer to hold chunks of incoming data

int const BUFFER_WORDS = 1;    // size of the buffer, in chars (2 bytes/char)
int const BUFFER_BYTES = 2*BUFFER_WORDS;
unsigned short buffer[BUFFER_WORDS];
