#ifndef TEXT_FILE_H_
#define TEXT_FILE_H_

#include <stdio.h>

class text_file
{
public:

	text_file(const char* fileName) : lineNumber(0)
	{
		file = fopen(fileName, "rt");
		if(file == 0)
			throw std::runtime_error("Failed to open file.");
	}

	~text_file()
	{
		fclose(file);
	}

	bool readLine(char *buffer)
	{
		if(fgets(buffer, readBufferSize, file)) {
			++lineNumber;
			return true;
		} else
			return false;
	}

	bool at_end() const
	{
		return feof(file) != 0;
	}

protected:
	static const size_t readBufferSize = 0x1000;
	FILE *file;
	size_t lineNumber;

};

#endif /* TEXT_FILE_H_ */
