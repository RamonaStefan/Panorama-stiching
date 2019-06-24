import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferInt;
import java.io.*;

import static java.lang.Float.MAX_VALUE;
import static java.lang.Float.MIN_VALUE;
import static java.lang.Math.ceil;
import static java.lang.Math.pow;
import static java.lang.StrictMath.*;


public class MyImage {
    public int width = 1;
    public int height = 1;
    public int noChannels;
    public float[][][] colorPattern;
    public BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_3BYTE_BGR);
    public BufferedImage finalImage;
    public static File filename;
    public static final float TWOPI = (float) 6.2831853;


    //resize image
    float nn_interpolate(float x, float y, int c) {
        return get_pixel(round(x), round(y), c);
    }

    MyImage nn_resize(int w, int h) {
        MyImage resizedImage = make_image(w, h, noChannels);
        float a_y = (float) width / w;
        float b_y = (float) (-0.5 + a_y * 0.5);
        float a_x = (float) height / h;
        float b_x = (float) (-0.5 + a_x * 0.5);
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                float x = a_x * i + b_x;
                float y = a_y * j + b_y;
                for (int c = 0; c < noChannels; c++) {
                    float value = nn_interpolate(x, y, c);
                    resizedImage.set_pixel(i, j, c, value);
                }
            }
        }
        return resizedImage;
    }

    float bilinear_interpolate(float x, float y, int c) {
        //floor x and y
        int xi = (int) floor(x);
        int yi = (int) floor(y);
        //get the values for all the corners
        float v1 = get_pixel(xi, yi, c);
        float v2 = get_pixel(xi, yi + 1, c);
        float v3 = get_pixel(xi + 1, yi, c);
        float v4 = get_pixel(xi + 1, yi + 1, c);
        //distances
        float d1 = x - xi;
        float d2 = 1 - d1;
        float d3 = y - yi;
        float d4 = 1 - d3;

        float q1 = v1 * d2 + v2 * d1;
        float q2 = v3 * d2 + v4 * d1;
        float q = q1 * d4 + q2 * d3;

        return q;
    }

    MyImage bilinear_resize(int w, int h) {
        MyImage resizedImage = make_image(w, h, noChannels);
        float a_y = (float) width / w;
        float b_y = (float) (-0.5 + a_y * 0.5);
        float a_x = (float) height / h;
        float b_x = (float) (-0.5 + a_x * 0.5);
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                float x = a_x * i + b_x;
                float y = a_y * j + b_y;
                for (int k = 0; k < noChannels; k++) {
                    float value = bilinear_interpolate(x, y, k);
                    resizedImage.set_pixel(i, j, k, value);
                }
            }
        }
        return resizedImage;
    }


    //filter image

    void l1_normalize() {
        float total = 0;
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                for (int k = 0; k < noChannels; k++) {
                    total += get_pixel(i, j, k);
                }
            }
        }
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                for (int k = 0; k < noChannels; k++) {
                    set_pixel(i, j, k, get_pixel(i, j, k) / total);
                }
            }
        }
    }


    MyImage make_box_filter(int w) {
        MyImage filterBox = make_image(w, w, 1);
        for (int i = 0; i < w; i++) {
            for (int j = 0; j < w; j++) {
                filterBox.set_pixel(i, j, 0, (float) 1 / (w * w));
            }
        }
        return filterBox;
    }

    MyImage convolve_image(MyImage filter, boolean preserve) {
        //filter should have one channel or the same number of channels as image
        assert (noChannels == filter.noChannels || filter.noChannels == 1);
        // filter should have an odd size
        assert (filter.height % 2 == 1);
        MyImage convolved_image;
        if (preserve == true) {
            convolved_image = make_image(width, height, noChannels);
        } else {
            convolved_image = make_image(width, height, 1);
        }
        for (int y = 0; y < convolved_image.height; y++) {
            for (int x = 0; x < convolved_image.width; x++) {
                if (preserve == true) {
                    for (int c = 0; c < noChannels; c++) {
                        int filter_c = c;
                        if (filter_c == 1) {
                            filter_c = 0;
                        }
                        float convolvedValue = apply_filter(filter, y, x, c, filter_c);
                        convolved_image.set_pixel(y, x, c, convolvedValue);
                    }
                } else {
                    float convolvedValue = 0;
                    for (int c = 0; c < noChannels; c++) {
                        int filter_c = c;
                        if (filter_c == 1) {
                            filter_c = 0;
                        }
                        convolvedValue += apply_filter(filter, y, x, c, filter_c);
                    }
                    convolved_image.set_pixel(y, x, 0, convolvedValue);
                }
            }
        }
        return convolved_image;

    }

    private float apply_filter(MyImage filter, int y, int x, int im_c, int filter_c) {
        float val = 0;
        int offset_x = filter.width / 2;
        int offset_y = filter.height / 2;
        for (int i = 0; i < filter.width; i++) {
            for (int j = 0; j < filter.height; j++) {
                float filter_pixel = filter.get_pixel(j, i, filter_c);
                float im_pixel = get_pixel(y - offset_y + j, x - offset_x + i, im_c);
                val += (filter_pixel * im_pixel);
            }
        }
        return val;
    }

    MyImage make_highpass_filter() {
        MyImage highpass = make_image(3, 3, 1);
        highpass.set_pixel(0, 1, 0, -1);
        highpass.set_pixel(1, 0, 0, -1);
        highpass.set_pixel(1, 1, 0, 4);
        highpass.set_pixel(1, 2, 0, -1);
        highpass.set_pixel(2, 1, 0, -1);
        return highpass;
    }

    MyImage make_sharpen_filter() {
        MyImage sharpen = make_image(3, 3, 1);
        sharpen.set_pixel(0, 1, 0, -1);
        sharpen.set_pixel(1, 0, 0, -1);
        sharpen.set_pixel(1, 1, 0, 5);
        sharpen.set_pixel(1, 2, 0, -1);
        sharpen.set_pixel(2, 1, 0, -1);
        return sharpen;
    }

    MyImage make_emboss_filter() {
        MyImage emboss = make_image(3, 3, 1);
        emboss.set_pixel(0, 0, 0, -2);
        emboss.set_pixel(0, 1, 0, -1);
        emboss.set_pixel(1, 0, 0, -1);
        emboss.set_pixel(1, 1, 0, 1);
        emboss.set_pixel(1, 2, 0, 1);
        emboss.set_pixel(2, 1, 0, 1);
        emboss.set_pixel(2, 2, 0, 2);
        return emboss;
    }

    MyImage make_gaussian_filter(float sigma) {
        double constant = 1 / (TWOPI * sigma * sigma);
        int sizeKernel = (int) ceil(6 * sigma);
        if (sizeKernel % 2 == 0) {
            sizeKernel++;
        }
        int offset = sizeKernel / 2;
        MyImage gauss = make_image(sizeKernel, sizeKernel, 1);
        for (int i = 0; i < sizeKernel; i++) {
            for (int j = 0; j < sizeKernel; j++) {
                double coords = -(pow(i - offset, 2) + pow(j - offset, 2));
                double val = constant * exp(coords / (double) (2 * sigma * sigma));
                gauss.set_pixel(i, j, 0, (float) val);
            }
        }

        gauss.l1_normalize();
        return gauss;
    }

    MyImage make_gx_filter() {

        MyImage gx = make_image(3, 3, 1);
        gx.set_pixel(0, 0, 0, -1);
        gx.set_pixel(0, 2, 0, 1);
        gx.set_pixel(1, 0, 0, -2);
        gx.set_pixel(1, 2, 0, 2);
        gx.set_pixel(2, 0, 0, -1);
        gx.set_pixel(2, 2, 0, 1);
        return gx;
    }

    MyImage make_gy_filter() {
        MyImage gy = make_image(3, 3, 1);
        gy.set_pixel(0, 0, 0, -1);
        gy.set_pixel(0, 1, 0, -2);
        gy.set_pixel(0, 2, 0, -1);
        gy.set_pixel(2, 0, 0, 1);
        gy.set_pixel(2, 1, 0, 2);
        gy.set_pixel(2, 2, 0, 1);
        return gy;
    }

    void feature_normalize() {
        float min = MAX_VALUE;
        float max = MIN_VALUE;
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                for (int k = 0; k < noChannels; k++) {
                    float pixel = get_pixel(i, j, k);
                    if (pixel < min) {
                        min = pixel;
                    }
                    if (pixel > max) {
                        max = pixel;
                    }
                }
            }
        }
        float range = max - min;
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                for (int c = 0; c < noChannels; c++) {
                    if (range != 0) {
                        float pixel = get_pixel(i, j, c);
                        pixel = (pixel - min) / range;
                        set_pixel(i, j, c, pixel);
                    } else {
                        set_pixel(i, j, c, 0);
                    }
                }

            }

        }
    }

    MyImage add_MyImage(MyImage b) {
        assert (noChannels == b.noChannels && height == b.height && width == b.width);
        MyImage result = make_image(width, height, noChannels);
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                for (int k = 0; k < noChannels; k++) {
                    float pixel = get_pixel(i, j, k);
                    pixel = pixel + b.get_pixel(i, j, k);
                    result.set_pixel(i, j, k, pixel);
                }
            }
        }
        return result;
    }

    MyImage sub_MyImage(MyImage b) {
        assert (noChannels == b.noChannels && height == b.height && width == b.width);
        MyImage result = make_image(width, height, noChannels);
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                for (int k = 0; k < noChannels; k++) {
                    float pixel = get_pixel(i, j, k);
                    pixel = pixel - b.get_pixel(i, j, k);
                    result.set_pixel(i, j, k, pixel);
                }
            }
        }
        return result;
    }

    //the gradient magnitude and direction
    MyImage[] sobel_image() {
        MyImage gx_filter = make_gx_filter();
        MyImage gy_filter = make_gy_filter();
        MyImage gx = convolve_image(gx_filter, false);
        MyImage gy = convolve_image(gy_filter, false);
        MyImage[] result = new MyImage[2];
        result[0] = make_image(width, height, 1);
        result[1] = make_image(width, height, 1);


        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                float pixel_gx = gx.get_pixel(i, j, 0);
                float pixel_gy = gx.get_pixel(i, j, 0);
                float val_mag = (float) sqrt(pow(pixel_gx, 2) + pow(pixel_gy, 2));
                float val_theta = (float) atan2(pixel_gy, pixel_gx);
                result[0].set_pixel(i, j, 0, val_mag);
                result[1].set_pixel(i, j, 0, val_theta);

            }
        }
        return result;
    }

    public MyImage(int noChannels, File filename) {

        this.filename = filename;
        try {
            image = ImageIO.read(filename);
        } catch (IOException e) {
            System.out.println("Could not find any photo with that name.");
        }
        width = image.getWidth();
        height = image.getHeight();
        this.noChannels = noChannels;
        colorPattern = new float[height][width][noChannels];
    }

    public MyImage() {
        this(3, new File("resources/dog.jpg"));
    }

    // Basic operations

    private float get_pixel(int r, int c, int ch) {
        if (r < 0) r = 0;
        if (r >= height) r = height - 1;
        if (c < 0) c = 0;
        if (c >= width) c = width - 1;
        if (ch < 0) ch = 0;
        if (ch >= noChannels) ch = noChannels - 1;
        return colorPattern[r][c][ch];
    }

    void set_pixel(int r, int c, int ch, float v) {
        if (c >= width || r >= height || ch >= noChannels || c < 0 || r < 0 || ch < 0) {
            return;
        }
        colorPattern[r][c][ch] = v;
    }


    MyImage copy_image() {
        MyImage im = make_image(width, height, noChannels);
        im.image = new BufferedImage(im.width, im.height, BufferedImage.TYPE_INT_RGB);
        im.finalImage = finalImage;
        im.filename = filename;
        for (int c = 0; c < noChannels; ++c) {
            for (int y = 0; y < height; ++y) {
                for (int x = 0; x < width; ++x) {
                    float pixel = get_pixel(y, x, c);
                    im.set_pixel(y, x, c, pixel);
                }
            }
        }
        return im;
    }

    // Loading and saving

    private MyImage make_image(int w, int h, int c) {
        MyImage im = new MyImage();
        im.height = h;
        im.width = w;
        im.noChannels = c;
        im.colorPattern = new float[h][w][c];
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                for (int k = 0; k < c; k++) {
                    im.set_pixel(i, j, k, 0);
                }
            }
        }
        return im;
    }

    void load_image() {
        // Create an empty BufferedImage object
        finalImage = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
        int[] temporary_array;
        float[] temporary_array_float = new float[width * height];

        // Draw Image into BufferedImage
        Graphics g = finalImage.getGraphics();
        g.drawImage(image, 0, 0, null);

        // Convert the BufferedImage to numeric pixel representation.
        DataBufferInt dataBufferInt = (DataBufferInt) finalImage.getRaster().getDataBuffer();
        temporary_array = dataBufferInt.getData();
        for (int i = 0; i < temporary_array.length; i++) {
            temporary_array_float[i] = (float) temporary_array[i];

        }

        colorPattern = convert13Dim(temporary_array_float, width, height, noChannels);
    }

    void save_image(String name) {
        String buff;
        buff = name + ".jpg";
        OutputStream fOut = null;
        try {
            fOut = new FileOutputStream(new File(buff));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        finalImage = new BufferedImage(width, height, BufferedImage.TYPE_3BYTE_BGR);

        for (int r = 0; r < height; r++) {
            for (int c = 0; c < width; c++) {
                if (noChannels > 1) {
                    int value = round(255 * colorPattern[r][c][0]) << 16 | round(255 * colorPattern[r][c][1]) << 8 | round(255 * colorPattern[r][c][2]);
                    finalImage.setRGB(c, r, value);
                } else {

                    int value = round(255 * colorPattern[r][c][0]) << 16 | round(255 * colorPattern[r][c][0]) << 8 | round(255 * colorPattern[r][c][0]);
                    finalImage.setRGB(c, r, value);
                }
            }
        }


        try {
            assert fOut != null;
            fOut.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

        try {
            // Write image
            ImageIO.write(finalImage, "jpg", new File(buff));
        } catch (IOException e) {
            System.out.println("Could not write the photo");
            System.err.println(e.getMessage());
        }
    }

    private float[][][] convert13Dim(float[] pixels1D, int Width, int Height, int noChannels) {
        float[][][] array3D = new float[Height][Width][noChannels];                  // New 3D array populated with color data

        for (int r = 0; r < Height; r++) {
            // Extract a row of pixel data into a temporary array of ints
            float[] temporaryArray = new float[Width];
            // Move the data into the 3D array.
            System.arraycopy(pixels1D, r * Width, temporaryArray, 0, Width);
            // Use the bitwise AND and bitwise right shift operations to mask all but the correct set of eight bits.
            for (int c = 0; c < Width; c++) {
                //Red data
                array3D[r][c][0] = (float) (((int) temporaryArray[c] >> 16) & 0xFF) / 255;
                //Green data
                array3D[r][c][1] = (float) (((int) temporaryArray[c] >> 8) & 0xFF) / 255;
                //Blue data
                array3D[r][c][2] = (float) (((int) temporaryArray[c]) & 0xFF) / 255;
            }
        }
        return array3D;
    }

    void shift_image(int c, float v) {
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                // get current value
                float val = get_pixel(y, x, c);
                // shift value by v
                val += v;
                // set pixel to new value
                set_pixel(y, x, c, val);
            }
        }
    }

    void clamp_image() {
        for (int c = 0; c < noChannels; ++c) {
            for (int y = 0; y < height; ++y) {
                for (int x = 0; x < width; ++x) {
                    float val = get_pixel(y, x, c);
                    if (val < 0) {
                        val = 0;
                    } else if (val > 1) {
                        val = 1;
                    }
                    set_pixel(y, x, c, val);
                }
            }
        }
    }

    void scale_image(int c, float v) {
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                // get current value
                float val = get_pixel(y, x, c);
                // shift value by v
                val *= v;
                // set pixel to new value
                set_pixel(y, x, c, val);
            }
        }
    }


    //panorama
    // Comparator for matches
    int match_compare(Match a, Match b) {
        if (a.distance < b.distance) return -1;
        else if (a.distance > b.distance) return 1;
        else return 0;
    }

    // Helper function to create 2d points.
    Point make_point(float x, float y) {
        Point p = new Point();
        p.x = x;
        p.y = y;
        return p;
    }

    // Place two images side by side on canvas, for drawing matching pixels.
    MyImage both_images(MyImage b) {
        MyImage both = make_image(width + b.width, height > b.height ? height : b.height, noChannels > b.noChannels ? noChannels : b.noChannels);
        for (int k = 0; k < noChannels; ++k) {
            for (int j = 0; j < height; ++j) {
                for (int i = 0; i < width; ++i) {
                    both.set_pixel(j, i, k, get_pixel(j, i, k));
                }
            }
        }
        for (int k = 0; k < b.noChannels; ++k) {
            for (int j = 0; j < b.height; ++j) {
                for (int i = 0; i < b.width; ++i) {
                    both.set_pixel(j, i + width, k, b.get_pixel(j, i, k));
                }
            }
        }
        return both;
    }


    // Find corners, match them, and draw them between two images.
    MyImage find_and_draw_matches(MyImage b, float sigma, float thresh, int nms) {
        int an = 0;
        int bn = 0;
        Descriptor[] ad = harris_corner_detector(sigma, thresh, nms);
        Descriptor[] bd = b.harris_corner_detector(sigma, thresh, nms);
        an = ad.length;
        bn = bd.length;
        Match[] m = match_descriptors(ad, bd);

        mark_corners(ad, an);
        b.mark_corners(bd, bn);
        MyImage lines = draw_matches(b, m, 0);

        return lines;
    }

    // Calculates L1 distance between to floating point arrays.
    float l1_distance(float[] a, float[] b) {
        float sum = 0;
        int an = a.length;
        int bn = b.length;
        assert (an != bn);
        int i;
        for (i = 0; i < bn; i++) {
            sum += abs(a[i] - b[i]);
        }
        return sum;
    }


    public void quickSort(Match[] arr, int low, int high) {
        if (arr == null || arr.length == 0)
            return;
        if (low >= high)
            return;
        // pick the pivot
        int middle = (low + high) / 2;
        Match pivot = arr[middle];
        // make left < pivot and right > pivot
        int i = low, j = high;
        while (i <= j) {
            while (match_compare(arr[i], pivot) == -1) {
                i++;
            }
            while (match_compare(arr[j], pivot) == 1) {
                j--;
            }
            if (i <= j) {
                Match temp = arr[i];
                arr[i] = arr[j];
                arr[j] = temp;
                i++;
                j--;
            }
        }
        // recursively sort two sub parts
        if (low < j)
            quickSort(arr, low, j);
        if (high > i)
            quickSort(arr, i, high);
    }

    // Finds best matches between descriptors of two images.
    Match[] match_descriptors(Descriptor[] a, Descriptor[] b) {
        // We will have at most an matches.
        int an = a.length;
        int bn = b.length;
        assert (an != bn);
        Match[] m = new Match[an];
        for (int i = 0; i < an; i++) {
            m[i] = new Match();
            m[i].made_zero();
        }
        float l1Distance;
        float maxDistance;
        for (int j = 0; j < an; ++j) {
            maxDistance = MAX_VALUE;
            for (int i = 0; i < bn; i++) {
                l1Distance = l1_distance(a[j].data, b[i].data);
                if (l1Distance < maxDistance) {
                    m[j].ai = j;
                    m[j].bi = i;
                    m[j].p = a[j].p;
                    m[j].q = b[i].p;
                    m[j].distance = l1Distance;
                    maxDistance = l1Distance;
                }
            }
        }

        int count = 0;
        int count2 = 0;
        int[] seen = new int[min(an, bn)];
        quickSort(m, 0, m.length - 1);
        boolean exists;
        for (int i = 0; i < m.length; i++) {
            exists = false;
            for (int j = 0; j < count; j++) {
                if (m[i].bi == seen[j]) {
                    exists = true;
                    break;
                }
            }
            if (exists == true) {
                Match aux;
                aux = m[i];
                for (int k = i; k < m.length - 1; k++) {
                    m[k] = m[k + 1];
                }
                m[m.length - 1] = aux;
                i--;
                count2++;
            } else {
                seen[count] = m[i].bi;
                count++;

            }
            if (count + count2 == m.length) {
                break;
            }
        }
        Match[] mFinal = new Match[count];
        for (int i = 0; i < count; i++) {
            mFinal[i] = m[i];
        }
        return mFinal;
    }

    // Apply a projective transformation to a point.
    Point project_point(Matrix H, Point p) {
        Matrix c = new Matrix(3, 1);
        c.data[0][0] = p.x;
        c.data[1][0] = p.y;
        c.data[2][0] = 1;
        Matrix homography = H.matrix_mult_matrix(c);
        float x = 0, y = 0;
        if (H.cols == 3 && H.rows == 3) {
            x = (float) (homography.data[0][0] / homography.data[2][0]);
            y = (float) (homography.data[1][0] / homography.data[2][0]);
        }
        Point q = make_point(x, y);
        return q;

    }

    // Calculate L2 distance between two points.
    float point_distance(Point p, Point q) {
        // TODO: should be a quick one.
        float error = (float) sqrt(pow(p.x - q.x, 2) + pow(p.y - q.y, 2));
        return error;
    }

    // Count number of inliers in a set of matches. Should also bring inliers
    int model_inliers(Matrix H, Match[] m, int n, float thresh) {
        int count = 0;
        int count2 = 0;
        for (int i = 0; i < n; i++) {
            Point pointX = project_point(H, m[i].p);
            Point pointY = m[i].q;
            float distance = point_distance(pointX, pointY);
            if (distance < thresh) {
                count++;
            } else {
                count2++;
                Match aux;
                aux = m[i];
                for (int k = i; k < n - 1; k++) {
                    m[k] = m[k + 1];
                }
                m[m.length - 1] = aux;
                i--;
            }
            if (count + count2 == n) {
                break;
            }
        }
        return count;
    }


//    float calculateMaxError(Matrix H, Match[] m, int n) {
//        float maxDistance = 0;
//        for (int i = 0; i < n; i++) {
//            Point pointX = project_point(H, m[i].p);
//            Point pointY = m[i].q;
//            float distance = point_distance(pointX, pointY);
//            maxDistance += distance;
//        }
//        return maxDistance / n;
//    }


    // Randomly shuffle matches for RANSAC.
    void randomize_matches(Match[] m, int n) {
        // TODO: implement Fisher-Yates to shuffle the array.
        for (int i = 0; i < n; i++) {
            int randomNumber = (int) (random() * i);
            Match aux;
            aux = m[i];
            m[i] = m[randomNumber];
            m[randomNumber] = aux;
        }

    }

    // Computes homography between two images given matching pixels.
    Matrix compute_homography(Match[] matches, int n) {
        Matrix M = new Matrix(n * 2, 8);
        Matrix b = new Matrix(n * 2, 1);

        for (int i = 0; i < n; ++i) {
            double x = matches[i].p.x;
            double xp = matches[i].q.x;
            double y = matches[i].p.y;
            double yp = matches[i].q.y;
            // TODO: fill in the matrices M and b.
            b.data[2 * i][0] = xp;
            b.data[2 * i + 1][0] = yp;
            M.data[2 * i][0] = x;
            M.data[2 * i][1] = y;
            M.data[2 * i][2] = 1;
            M.data[2 * i][6] = -x * xp;
            M.data[2 * i][7] = -y * xp;
            M.data[2 * i + 1][3] = x;
            M.data[2 * i + 1][4] = y;
            M.data[2 * i + 1][5] = 1;
            M.data[2 * i + 1][6] = -x * yp;
            M.data[2 * i + 1][7] = -y * yp;

        }

        Matrix a = M.solve_system(b);
        // If a solution can't be found, return empty matrix;
        Matrix none = new Matrix();
        if (a.data == null) return none;
        Matrix H = new Matrix(3, 3);

        if (a.cols == 1 && a.rows == 8) {
            //line 1
            H.data[0][0] = a.data[0][0];
            H.data[0][1] = a.data[1][0];
            H.data[0][2] = a.data[2][0];

            //line 2
            H.data[1][0] = a.data[3][0];
            H.data[1][1] = a.data[4][0];
            H.data[1][2] = a.data[5][0];

            //line 3
            H.data[2][0] = a.data[6][0];
            H.data[2][1] = a.data[7][0];
            H.data[2][2] = 1;
        }
        return H;
    }

    // Perform RANdom SAmple Consensus to calculate homography for noisy matches.
    Matrix RANSAC(Match[] m, int n, float thresh, int k, int cutoff) {
        int bestfit = 0;
        Matrix aux = new Matrix();
        Matrix Hb = aux.make_translation_homography(256, 0);
        for (int i = 0; i < k; i++) {
            randomize_matches(m, m.length);
            Matrix H = compute_homography(m, n);
            int inliers = model_inliers(H, m, m.length, thresh);
            if (inliers > bestfit) {
                H = compute_homography(m, inliers);
                bestfit = inliers;
                Hb = H;
                if (inliers > cutoff) {
                    return Hb;
                }
            }
        }
        return Hb;
    }

    // Stitches two images together using a projective transformation.
    MyImage combine_images(MyImage b, Matrix H) {
        Matrix Hinv = H.matrix_invert();
        Matrix none = new Matrix();
        // Project the corners of image b into image a coordinates.
        Point c1 = b.project_point(Hinv, make_point(0, 0));
        Point c2 = b.project_point(Hinv, make_point(b.width - 1, 0));
        Point c3 = b.project_point(Hinv, make_point(0, b.height - 1));
        Point c4 = b.project_point(Hinv, make_point(b.width - 1, b.height - 1));

        // Find top left and bottom right corners of image b warped into image a.
        Point topleft = new Point();
        Point botright = new Point();
        botright.x = max(c1.x, max(c2.x, max(c3.x, c4.x)));
        botright.y = max(c1.y, max(c2.y, max(c3.y, c4.y)));
        topleft.x = min(c1.x, min(c2.x, min(c3.x, c4.x)));
        topleft.y = min(c1.y, min(c2.y, min(c3.y, c4.y)));

        // Find how big our new image should be and the offsets from image a.
        int dx = (int) min(0, topleft.x);
        int dy = (int) min(0, topleft.y);
        int w = (int) (max(width, botright.x) - dx);
        int h = (int) (max(height, botright.y) - dy);

        MyImage c = make_image(w, h, noChannels);
        // Paste image a into the new image offset by dx and dy.
      /*  for (int k = 0; k < noChannels; ++k) {
            for (int j = 0; j < height; ++j) {
                for (int i = 0; i < width; ++i) {
                    // TODO: fill in.
                    c.set_pixel(j - dy, i - dx, k, get_pixel(j, i, k));
                }
            }
        }
        for (int i = (int) topleft.x; i <= botright.x; i++) {
            for (int j = (int) topleft.y; j <= botright.y; j++) {
                for (int k = 0; k < c.noChannels; k++) {
                    Point p = project_point(H, new Point(i, j));
                    if (p.y >= 0 && p.y < b.height && p.x >= 0 && p.x < b.width) {
                        float value = b.bilinear_interpolate(p.y, p.x, k);
                        c.set_pixel(j - dy, i - dx, k, value);
                    }
                }
            }
        }
     */
        for (int k = 0; k < c.noChannels; ++k) {
            for (int j = 0; j < c.height; ++j) {
                for (int i = 0; i < c.width; ++i) {
                    if (j + dy >= 0 && j + dy < height && i + dx >= 0 && i + dx < width) {
                        c.set_pixel(j, i, k, get_pixel(j + dy, i + dx, k));
                    }
                    if (i + dx >= topleft.x && i + dx <= botright.x && j + dy >= 0 && j + dy <= botright.y) {
                        Point p = project_point(H, new Point(i + dx, j + dy));
                        if (p.y >= 0 && p.y < b.height && p.x >= 0 && p.x < b.width) {
                            float value = b.bilinear_interpolate(p.y, p.x, k);
                            c.set_pixel(j, i, k, value);
                        }
                    }
                }
            }
        }

        return c;
    }

    // Create a panoramam between two images.
    MyImage panorama_image(MyImage b, float sigma, float thresh, int nms, float inlier_thresh, int iters, int cutoff) {
        int an = 0;
        int bn = 0;

        // Calculate corners and descriptors
        Descriptor[] ad = harris_corner_detector(sigma, thresh, nms);
        Descriptor[] bd = b.harris_corner_detector(sigma, thresh, nms);
        an = ad.length;
        bn = bd.length;
        if (an == 0 || bn == 0) {
            return this;
        }

        // Find matches
        Match[] m = match_descriptors(ad, bd);
        if (m.length == 0) {
            return this;
        }

        // Run RANSAC to find the homography
        Matrix H;
        Matrix none = new Matrix();
        do {
            H = RANSAC(m, 4, inlier_thresh, iters, cutoff);
        } while (H.matrix_invert() == none);

        if (true) {
            // Mark corners and matches between images
            //mark_corners(ad, an);
            // b.mark_corners(bd, bn);
            MyImage inlier_matches = draw_inliers(b, H, m, m.length, inlier_thresh);
            inlier_matches.save_image("out/inliers");
        }

        // Stitch the images together with the homography
        MyImage comb = combine_images(b, H);
        return comb;
    }

    // Project an image onto a cylinder.
    MyImage cylindrical_project(float f) {
        //TODO: project image onto a cylinder
        MyImage c = make_image(width, height, noChannels);
        int xc = width / 2;
        int yc = height / 2;
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                for (int k = 0; k < noChannels; k++) {
                    float theta = (i - xc) / f;
                    float h = (j - yc) / f;
                    float X = (float) sin(theta);
                    float Y = h;
                    float Z = (float) cos(theta);
                    float x = f * X / Z + xc;
                    float y = f * Y / Z + yc;
                    //   if (x >= 0 && x < width && y >=0 && y < height) {
                    float value = bilinear_interpolate(y, x, k);
                    c.set_pixel(j, i, k, value);
                    //  }
                }
            }
        }

        return c;
    }


    //harris
    // Calculate the structure matrix of an image.
    MyImage structure_matrix(float sigma) {
        MyImage gray = rgb_to_grayscale();
        for (int i = 0; i < gray.height; i++) {
            for (int j = 0; j < gray.width; j++) {
                for (int k = 0; k < gray.noChannels; k++) {
                    gray.colorPattern[i][j][k] *= 255;
                }
            }
        }
        MyImage gauss = make_gaussian_filter(sigma);
        MyImage gx_filter = make_gx_filter();
        MyImage gy_filter = make_gy_filter();
        MyImage Gx = gauss.convolve_image(gx_filter, false);
        MyImage Gy = gauss.convolve_image(gy_filter, false);
        MyImage Ix = gray.convolve_image(Gx, false);
        MyImage Iy = gray.convolve_image(Gy, false);

        MyImage S = make_image(width, height, 3);
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                S.set_pixel(i, j, 0, Ix.get_pixel(i, j, 0) * Ix.get_pixel(i, j, 0));
                S.set_pixel(i, j, 1, Iy.get_pixel(i, j, 0) * Iy.get_pixel(i, j, 0));
                S.set_pixel(i, j, 2, Ix.get_pixel(i, j, 0) * Iy.get_pixel(i, j, 0));
            }
        }


        S = S.convolve_image(gauss, true);
        return S;
    }

    // Estimate the cornerness of each pixel given a structure matrix S.
    MyImage cornerness_response() {
        MyImage R = make_image(width, height, 1);
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                double det = (get_pixel(i, j, 0) * get_pixel(i, j, 1)) - (get_pixel(i, j, 2) * get_pixel(i, j, 2));
                double trace = get_pixel(i, j, 0) + get_pixel(i, j, 1);
                R.set_pixel(i, j, 0, (float) ((det - 0.06 * trace * trace)));

            }
        }
        return R;
    }

    // Perform non-max supression on an image of feature responses.
    MyImage nms_image(int w) {
        MyImage r = copy_image();
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                for (int k = 0; k < noChannels; k++) {
                    check(i, j, k, w, r);
                }
            }
        }
        return r;
    }

    void check(int x, int y, int k, int w, MyImage r) {
        for (int i = -w; i <= w; i++) {
            for (int j = -w; j <= w; j++) {
                if (get_pixel(x + i, y + j, k) > r.get_pixel(x, y, k)) {
                    r.set_pixel(x, y, k, -999999);
                }
            }
        }
    }

    // Perform harris corner detection and extract features from the corners.
    Descriptor[] harris_corner_detector(float sigma, float thresh, int nms) {
        // Calculate structure matrix
        MyImage S = structure_matrix(sigma);

        // Estimate cornerness
        MyImage R = S.cornerness_response();

        // Run NMS on the responses
        MyImage Rnms = R.nms_image(nms);

        int count = 0;
        for (int i = 0; i < Rnms.height; i++) {
            for (int j = 0; j < Rnms.width; j++) {
                for (int k = 0; k < Rnms.noChannels; k++) {
                    if (Rnms.get_pixel(i, j, k) / 255 / 255 > thresh) {
                        count++;
                    }
                }
            }
        }
        int index = 0;
        Descriptor[] d = new Descriptor[count];
        for (int i = 0; i < Rnms.height; i++) {
            for (int j = 0; j < Rnms.width; j++) {
                for (int k = 0; k < Rnms.noChannels; k++) {
                    if (Rnms.get_pixel(i, j, k) / 255 / 255 > thresh) {///pow(10,3)) {
                        d[index++] = describe_index(i * Rnms.width + j);
                    }
                }
            }
        }
        return d;
    }

    // Create a feature descriptor for an index in an image.
    Descriptor describe_index(int i) {
        int w = 5;
        Descriptor d = new Descriptor();
        d.p.x = i % width;
        d.p.y = i / width;
        d.data = new float[w * w * noChannels];
        for (int l = 0; l < w * w * noChannels; l++) {
            d.data[l] = 0;
        }
        d.n = w * w * noChannels;
        int count = 0;
        // If you want you can experiment with other descriptors
        // This subtracts the central value from neighbors
        // to compensate some for exposure/lighting changes.
        for (int c = 0; c < noChannels; ++c) {
            float cval = this.get_pixel(i / width, i % width, c);
            for (int dx = -w / 2; dx < (w + 1) / 2; ++dx) {
                for (int dy = -w / 2; dy < (w + 1) / 2; ++dy) {
                    float val = this.get_pixel(i / width + dy, i % width + dx, c);
                    d.data[count++] = cval - val;
                }
            }
        }
        return d;
    }

    // Find and draw corners on an image.
    void detect_and_draw_corners(float sigma, float thresh, int nms) {
        Descriptor[] d = harris_corner_detector(sigma, thresh, nms);
        mark_corners(d, d.length);
    }

    // Marks the spot of a point in an image.
    void mark_spot(Point p) {
        int x = (int) p.x;
        int y = (int) p.y;
        int i;
        for (i = -9; i < 10; ++i) {
            set_pixel(y, x + i, 0, 1);
            set_pixel(y + i, x, 0, 1);
            set_pixel(y, x + i, 1, 0);
            set_pixel(y + i, x, 1, 0);
            set_pixel(y, x + i, 2, 1);
            set_pixel(y + i, x, 2, 1);
        }
    }

    // Marks corners denoted by an array of descriptors.
    void mark_corners(Descriptor[] d, int n) {
        int i;
        for (i = 0; i < n; ++i) {
            mark_spot(d[i].p);
        }
    }

    // Draws lines between matching pixels in two images.
    MyImage draw_matches(MyImage b, Match[] matches, int inliers) {
        MyImage both = both_images(b);
        for (int i = 0; i < matches.length; ++i) {
            int bx = (int) matches[i].p.x;
            int ex = (int) matches[i].q.x;
            int by = (int) matches[i].p.y;
            int ey = (int) matches[i].q.y;
            for (int j = bx; j < ex + width; ++j) {
                int r = (int) ((float) (j - bx) / (ex + width - bx) * (ey - by) + by);
                both.set_pixel(r, j, 0, i < inliers ? 0 : 1);
                both.set_pixel(r, j, 1, i < inliers ? 1 : 0);
                both.set_pixel(r, j, 2, 0);
            }
        }
        return both;
    }

    // Draw the matches with inliers in green between two images.
    MyImage draw_inliers(MyImage b, Matrix H, Match[] m, int n, float thresh) {
        int inliers = model_inliers(H, m, n, thresh);
        MyImage lines = draw_matches(b, m, inliers);
        return lines;
    }

    MyImage rgb_to_grayscale() {
        assert (noChannels == 3);
        MyImage gray = make_image(width, height, 1);
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                // get all three colors
                float r = get_pixel(y, x, 0);
                float g = get_pixel(y, x, 1);
                float b = get_pixel(y, x, 2);

                // calculate new value
                float val = (float) (0.299 * r + 0.587 * g + 0.114 * b);
                gray.set_pixel(y, x, 0, val);
            }
        }


        return gray;
    }
}
