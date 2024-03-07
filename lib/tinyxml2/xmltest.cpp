

#include "tinyxml2.h"
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <string>
#include <iostream>




void testFunction(){
	tinyxml2::XMLDocument xml_doc;

	tinyxml2::XMLError eResult = xml_doc.LoadFile("resources/myXML.xml");
	tinyxml2::XMLElement* root = xml_doc.FirstChildElement("IO");

	tinyxml2::XMLElement* inputMeshFileElt = root->FirstChildElement("InputMeshFile");
	tinyxml2::XMLElement* outputFilePathElt = root->FirstChildElement("OutputFilePath");

	std::string inputMeshFile(inputMeshFileElt->GetText()); 
	std::string outputFilePath(outputFilePathElt->GetText()); 

	std::cout << inputMeshFile << std::endl;
	std::cout << outputFilePath << std::endl;
	
}



int main( int argc, const char ** argv )
{
	testFunction();
	
	
	return 0;
}