#include <iostream>
#include <string>
#include <sstream>
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>

#include "fluid_dyn/fluid.h"
#include "utils/utils.h"
#include <conio.h>

using namespace std;

class controller
{
private:
    sf::Image image;
    sf::Image bg;
    sf::Image alpha_channel;
    sf::Texture texture;
    int x,y;
    float** initial_scene;
    fluid_space* space;

    void set_initial_scene()
    {
        initial_scene=Utils::allocate_matrix<float>(y,x);
        int i,j;
            for(i=0;i<y;i++)
                for(j=0;j<x;j++)
                {
                    initial_scene[i][j]=(float)(image.getPixel(j,i).r+image.getPixel(j,i).g+image.getPixel(j,i).b)/255/3;
                    //image.setPixel(j,i,sf::Color(initial_scene[i][j]*255,initial_scene[i][j]*255,initial_scene[i][j]*255));
                }
    }

    sf::Image get_Image_from_field(float** matrix)
    {
        sf::Image img;
        img.create(x,y);
        int i,j;
        sf::Color color;
        for(i=0;i<x;i++)
            for(j=0;j<y;j++)
                if(matrix[j][i]<0)
                    {
                        color=bg.getPixel(i,j);
                        color=sf::Color::White;
                        img.setPixel(i,j,color);
                    }
                else
                    if(matrix[j][i]>1)
                        img.setPixel(i,j,bg.getPixel(i,j));
                    else
                    {
                        color=bg.getPixel(i,j);

                        if(color.r+(1-matrix[j][i])*(alpha_channel.getPixel(i,j).r)>255)
                            color.r=255;
                        else color.r+=(1-matrix[j][i])*(alpha_channel.getPixel(i,j).r);

                        if(color.g+(1-matrix[j][i])*(alpha_channel.getPixel(i,j).g)>255)
                            color.g=255;
                        else color.g+=(1-matrix[j][i])*(alpha_channel.getPixel(i,j).g);

                        if(color.b+(1-matrix[j][i])*(alpha_channel.getPixel(i,j).b)>255)
                            color.b=255;
                        else color.b+=(1-matrix[j][i])*(alpha_channel.getPixel(i,j).b);

                        img.setPixel(i,j,color);
                    }
        return img;
    }



public:
    controller(char* bgf, char* filename, char* alpha)
    {
        if (!image.loadFromFile(filename))
            throw 1;
        if (!bg.loadFromFile(bgf))
            throw 1;
        if (!alpha_channel.loadFromFile(alpha))
            throw 1;

        x=image.getSize().x;
        y=image.getSize().y;
        space=new fluid_space(x,0.01,0.05,0.0000001);
        set_initial_scene();

        texture.loadFromImage(image);
    }

    void render()
    {

        sf::RenderWindow window(sf::VideoMode(x, y), "SFML works!");

        space->add_density(initial_scene);

        int frame=100;

        while (window.isOpen())
        {
            sf::Event event;
            while (window.pollEvent(event))
            {
                if (event.type == sf::Event::Closed)
                    window.close();
            }

            {
                window.clear();
                sf::Sprite sprite;
                float** density=space->get_density();
                texture.loadFromImage(get_Image_from_field(density));
                Utils::deallocate_matrix<float>(density,x);
                sprite.setTexture(texture);
                string path="C:\\Users\\Alex\\Desktop\\CDROM_GDC03\\bufny frames\\";
                stringstream ss;
                ss << frame;
                path.append(ss.str());
                path.append(".bmp");
                cout<< path;
                frame++;
                texture.copyToImage().saveToFile(path.c_str());
                sprite.setPosition(0,0);
                window.draw(sprite);
                window.display();
                sf::sleep(sf::milliseconds(00));
                space->step();
            }
        }
    }
    ~controller()
    {
        Utils::deallocate_matrix(initial_scene,x);
        delete space;
    }

};

