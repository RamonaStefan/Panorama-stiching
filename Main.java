import java.io.File;

public class Main {

    public static void main(String[] args) {
        MyImage object = new MyImage();
        object.load_image();

        // nearest resize
        //object = object.nn_resize(object.width*4, object.height*4);

        // bilinear resize
        //object = object.bilinear_resize(object.width*4, object.height*4);
       // object.save_image("out/final2");

        // box filter + resize
        //MyImage filter = object.make_box_filter(7);
        //MyImage convolvedImage = object.convolve_image(filter, true);
        //convolvedImage = convolvedImage.nn_resize(convolvedImage.width/7, convolvedImage.height/7);
        //convolvedImage.save_image("out/final4");

        // apply gaussian filter
        //MyImage gauss = object.make_gaussian_filter(2);
       //MyImage convolvedImage = object.convolve_image(gauss, true);
       //convolvedImage.save_image("out/final6");

        // hybrid images - low and high frequencies
        //MyImage gauss = object.make_gaussian_filter(2);
        //MyImage lfreq = object.convolve_image(gauss, true);
        //MyImage hfreq = object.sub_MyImage(lfreq);
        //MyImage reconstruct = lfreq.add_MyImage(hfreq);
        //lfreq.save_image("out/low");
        //hfreq.save_image("out/high");
        //reconstruct.save_image("out/reconstruct");

        //sobel
        //MyImage[] sobel = object.sobel_image();
        //MyImage mag = sobel[0];
        //mag.feature_normalize();
        //mag.save_image("out/final7");

        //corner detection
       // MyImage object2 = new MyImage(3, new File("resources/Rainier1.jpg"));
       //object2.load_image();
       //object2.detect_and_draw_corners(2,(float)500,3);
       //object2.save_image("out/corners");

//        //matching
//       MyImage object6 = new MyImage(3, new File("resources/Rainier1.jpg"));
//       MyImage object7 = new MyImage(3, new File("resources/Rainier2.jpg"));
//       object6.load_image();
//       object7.load_image();
//       MyImage m = object7.find_and_draw_matches(object6,2, 10, 3);
//       m.save_image("out/matches");

//        //panorama
//        // 2 photos
//        MyImage object2 = new MyImage(3, new File("resources/Rainier1.jpg"));
//        MyImage object3 = new MyImage(3, new File("resources/Rainier2.jpg"));
//       object2.load_image();
//       object3.load_image();
//       object2= object2.cylindrical_project(8000);
//       object3 = object3.cylindrical_project(8000);
//       MyImage panorama = object2.panorama_image(object3, 2, 10, 3, 10, 10000, 50);
//       panorama = panorama.cylindrical_project(8000);
//       //3 photos
//       MyImage object4 = new MyImage(3, new File("resources/Rainier5.jpg"));
//       object4.load_image();
//       object4 = object4.cylindrical_project(8000);
//       panorama = panorama.panorama_image(object4, 2, 10, 3, 10, 10000, 50);
//        panorama = panorama.cylindrical_project(8000);
//        // 4 photos
//        MyImage object5 = new MyImage(3, new File("resources/Rainier6.jpg"));
//        object5.load_image();
//        object5 = object5.cylindrical_project(8000);
//        panorama = panorama.panorama_image(object5, 2, 10, 3, 10, 10000, 50);
//       panorama.save_image("out/panorama");

       //5 photos
//        MyImage object6 = new MyImage(3, new File("resources/Rainier3.jpg"));
//        object6.load_image();
//        object6 = object6.cylindrical_project(8000);
//        panorama = panorama.panorama_image(object6, 2, 10, 3, 10, 10000, 50);
//        panorama.save_image("out/panorama");
//
//        //6 photos
//        MyImage object7 = new MyImage(3, new File("resources/Rainier3.jpg"));
//        object7.load_image();
//        object7 = object7.cylindrical_project(8000);
//        panorama = panorama.panorama_image(object7, 2, 10, 3, 10, 10000, 50);
//        panorama.save_image("out/panorama");




        //panorama
        // 2 photos
        MyImage object2 = new MyImage(3, new File("resources/field2.jpg"));
        MyImage object3 = new MyImage(3, new File("resources/field3.jpg"));
        object2.load_image();
        object3.load_image();
        object2= object2.cylindrical_project(6000);
        object3 = object3.cylindrical_project(6000);
        MyImage panorama = object2.panorama_image(object3, 2, 40, 3, 20, 10000, 50);
        panorama = panorama.cylindrical_project(6000);
        panorama.save_image("out/panorama");


        //panorama
        // 2 photos
//        MyImage object2 = new MyImage(3, new File("resources/sun3.jpg"));
//        MyImage object3 = new MyImage(3, new File("resources/sun4.jpg"));
//        object2.load_image();
//        object3.load_image();
//        object2= object2.cylindrical_project(6000);
//        object3 = object3.cylindrical_project(6000);
//        MyImage panorama = object2.panorama_image(object3, 2, 20, 3, 20, 10000, 50);
//        panorama = panorama.cylindrical_project(6000);
//       panorama.save_image("out/panorama");





    }
}
