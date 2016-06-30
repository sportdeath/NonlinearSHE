#include <CImg.h>

int main() {
  cimg_library::CImg<unsigned char> img("image.jpg");
  //img.fill(0);
  //unsigned char purple[] = {255,0,255};
  //img.draw_text(0,0,"Hello World", purple);
  img.display("wowee");
  return 0;
}
