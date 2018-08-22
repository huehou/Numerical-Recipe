#include <stdio.h>

/*  For Question 1(b)(i)
void fun(char *s[])
{
    for(int i = 0; i < 5; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            printf("%c", s[i][j]);
        }
        printf("\n");
    }
}

int main()
{
    char *s[5] = {"this", "that", "we\0\0", "!\0\0\0","guys"};
    fun(s);

    printf("\nPress ENTER to exit...\n");
    getchar();
    return 0;
} */

void fun(char s[][5])
{
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            printf("%c", s[i][j]);
        }
        printf("\n");
    }
}

int main()
{
    char t[4][5] = {"this", "that", "we", "!"};
    fun(t);

    printf("\nPress ENTER to exit...\n");
    getchar();
    return 0;
}