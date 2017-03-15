#include <iostream>

#include "src/controller.h"

using namespace std;

LRESULT CALLBACK WndProc(HWND, UINT, WPARAM, LPARAM);
char *OpenDialog(HWND);
HWND ghwndEdit;

void OpenDialog(HWND hwnd, char* filepath, char* filter)
{
    OPENFILENAME ofn;
    TCHAR szFile[MAX_PATH];

    ZeroMemory(&ofn, sizeof(ofn));
    ofn.lStructSize = sizeof(ofn);
    ofn.lpstrFile = szFile;
    ofn.lpstrFile[0] = '\0';
    ofn.hwndOwner = hwnd;
    ofn.nMaxFile = sizeof(szFile);
    ofn.lpstrFilter = TEXT(filter);
    ofn.nFilterIndex = 1;
    ofn.lpstrInitialDir = NULL;
    ofn.lpstrFileTitle = NULL;
    ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

    if(GetOpenFileName(&ofn))
    {
        //LoadFile(ofn.lpstrFile);
        CloseWindow(hwnd);
        //std::cout<<ofn.lpstrFile;
        strcpy(filepath,ofn.lpstrFile);
    }
    else
        filepath[0]=0;
}

int main()
{
    char filename[256],bgf[256],alpha[256];

    OpenDialog(ghwndEdit,filename,"Image Files\0""*.*\0");
    if(!filename[0])
        return 1;

    OpenDialog(ghwndEdit,bgf,"Image Files\0""*.*\0");
    if(!bgf[0])
        return 1;

    OpenDialog(ghwndEdit,alpha,"Image Files\0""*.*\0");
    if(!alpha[0])
        return 1;

    controller ctr=controller(bgf,filename,alpha);
    ctr.render();
    return 0;
}
