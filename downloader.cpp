#include <iostream>
#include <sstream>
#include <cmath>
#include <fstream>
#include <vector>
#include <thread>
#include <mutex>
#include <curl/curl.h>
#include <gdal_priv.h>
#include <cpl_conv.h>

class Downloader {
public:
    Downloader(int index, int count, const std::vector<std::string>& urls, std::vector<std::vector<char>>& datas)
        : index(index), count(count), urls(urls), datas(datas) {
    }

    ~Downloader() {
    }

    void download(const std::string& url) {
        CURL* curl = curl_easy_init();
        if (curl) {
            printf("%s\n", url.data());
            curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
            curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, &Downloader::WriteCallback);
            curl_easy_setopt(curl, CURLOPT_WRITEDATA, this);
            curl_easy_setopt(curl, CURLOPT_USERAGENT, "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/88.0.4324.150 Safari/537.36 Edg/88.0.705.68");

            CURLcode res = curl_easy_perform(curl);
            if (res != CURLE_OK) {
                std::cerr << "curl_easy_perform() failed: " << curl_easy_strerror(res) << std::endl;
            }
            curl_easy_cleanup(curl);
        }
    }

    static size_t WriteCallback(void* contents, size_t size, size_t nmemb, void* userp) {
        Downloader* downloader = static_cast<Downloader*>(userp);
        size_t totalSize = size * nmemb;
        printf("totalSize:%d\n", totalSize);
        downloader->datas[downloader->index].insert(downloader->datas[downloader->index].end(),
                                                    static_cast<char*>(contents),
                                                    static_cast<char*>(contents) + totalSize);
        return totalSize;
    }

    void operator()() {
        for (size_t i = index; i < urls.size(); i += count) {
            download(urls[i]);
        }
    }

private:
    int index;
    int count;
    const std::vector<std::string>& urls;
    std::vector<std::vector<char>>& datas;
};

// WGS-84 to Web Mercator
std::pair<double, double> wgsToMercator(double x, double y) {
    y = (y > 85.0511287798) ? 85.0511287798 : y;
    y = (y < -85.0511287798) ? -85.0511287798 : y;

    double x2 = x * 20037508.34 / 180;
    double y2 = log(tan((90 + y) * M_PI / 360)) / (M_PI / 180);
    y2 = y2 * 20037508.34 / 180;
    return std::make_pair(x2, y2);
}

// Web Mercator to WGS-84
std::pair<double, double> mercatorToWgs(double x, double y) {
    double x2 = x / 20037508.34 * 180;
    double y2 = y / 20037508.34 * 180;
    y2 = 180 / M_PI * (2 * atan(exp(y2 * M_PI / 180)) - M_PI / 2);
    return std::make_pair(x2, y2);
}

// Functions for coordinate conversion between GCJ-02 and WGS-84 (not fully translated)
// ...

// Get tile coordinates in Google Maps based on latitude and longitude of WGS-84
std::pair<int, int> wgsToTile(double j, double w, int z) {
    // auto isNum = [](double x) { return std::is_floating_point<double>(x); };

    // if (!isNum(j) || !isNum(w)) {
    //     throw std::invalid_argument("j and w must be int or float!");
    // }

    if ( z < 0 || z > 22) {
        throw std::invalid_argument("z must be int and between 0 to 22.");
    }

    if (j < 0) {
        j = 180 + j;
    } else {
        j += 180;
    }
    j /= 360;

    w = (w > 85.0511287798) ? 85.0511287798 : w;
    w = (w < -85.0511287798) ? -85.0511287798 : w;
    w = log(tan((90 + w) * M_PI / 360)) / (M_PI / 180);
    w /= 180;
    w = 1 - (w + 1) / 2;

    int num = std::pow(2, z);
    int x = std::floor(j * num);
    int y = std::floor(w * num);
    return std::make_pair(x, y);
}

// 在GCJ-02坐标和WGS-84坐标之间进行转换
// 中国大陆的所有公共地理数据需要使用GCJ-02坐标进行加密，引入了随机偏移
// 以下部分代码用于消除这种偏移
double transformLat(double x, double y) {
    double ret = -100.0 + 2.0 * x + 3.0 * y + 0.2 * y * y + 0.1 * x * y + 0.2 * sqrt(abs(x));
    ret += (20.0 * sin(6.0 * x * M_PI) + 20.0 * sin(2.0 * x * M_PI)) * 2.0 / 3.0;
    ret += (20.0 * sin(y * M_PI) + 40.0 * sin(y / 3.0 * M_PI)) * 2.0 / 3.0;
    ret += (160.0 * sin(y / 12.0 * M_PI) + 320 * sin(y * M_PI / 30.0)) * 2.0 / 3.0;
    return ret;
}

double transformLon(double x, double y) {
    double ret = 300.0 + x + 2.0 * y + 0.1 * x * x + 0.1 * x * y + 0.1 * sqrt(abs(x));
    ret += (20.0 * sin(6.0 * x * M_PI) + 20.0 * sin(2.0 * x * M_PI)) * 2.0 / 3.0;
    ret += (20.0 * sin(x * M_PI) + 40.0 * sin(x / 3.0 * M_PI)) * 2.0 / 3.0;
    ret += (150.0 * sin(x / 12.0 * M_PI) + 300.0 * sin(x / 30.0 * M_PI)) * 2.0 / 3.0;
    return ret;
}

// 计算坐标偏移
std::map<std::string, double> delta(double lat, double lon) {
    // Krasovsky 1940
    // a = 6378245.0, 1/f = 298.3
    // b = a * (1 - f)
    // ee = (a^2 - b^2) / a^2
    double a = 6378245.0; // a: 卫星椭球体坐标映射到平面地图坐标系的投影因子
    double ee = 0.00669342162296594323; // ee: 椭球体离心率
    double dLat = transformLat(lon - 105.0, lat - 35.0);
    double dLon = transformLon(lon - 105.0, lat - 35.0);
    double radLat = lat / 180.0 * M_PI;
    double magic = sin(radLat);
    magic = 1 - ee * magic * magic;
    double sqrtMagic = sqrt(magic);
    dLat = (dLat * 180.0) / ((a * (1 - ee)) / (magic * sqrtMagic) * M_PI);
    dLon = (dLon * 180.0) / (a / sqrtMagic * cos(radLat) * M_PI);
    
    std::map<std::string, double> result;
    result["lat"] = dLat;
    result["lon"] = dLon;
    return result;
}

// 判断是否在中国境内
bool outOfChina(double lat, double lon) {
    return (lon < 72.004 || lon > 137.8347) || (lat < 0.8293 || lat > 55.8271);
}

// 从GCJ-02坐标转换为WGS-84坐标
std::pair<double, double> gcj_to_wgs(double gcjLon, double gcjLat) {
    if (outOfChina(gcjLat, gcjLon)) {
        return std::make_pair(gcjLon, gcjLat);
    }
    std::map<std::string, double> d = delta(gcjLat, gcjLon);
    return std::make_pair(gcjLon - d["lon"], gcjLat - d["lat"]);
}

// 从WGS-84坐标转换为GCJ-02坐标
std::pair<double, double> wgs_to_gcj(double wgsLon, double wgsLat) {
    if (outOfChina(wgsLat, wgsLon)) {
        return std::make_pair(wgsLon, wgsLat);
    }
    std::map<std::string, double> d = delta(wgsLat, wgsLon);
    return std::make_pair(wgsLon + d["lon"], wgsLat + d["lat"]);
}

// Function to convert pixel coordinates to Mercator coordinates
std::map<std::string, std::pair<double, double>> pixlsToMercator(const std::map<std::string, std::pair<int, int>>& zb,int z) {
    int inx = zb.at("LT").first;  // left top
    int iny = zb.at("LT").second;
    int inx2 = zb.at("RB").first;  // right bottom
    int iny2 = zb.at("RB").second;
    double length = 20037508.3427892;
    int sum = 1 << z; // Equivalent to 2 raised to the power of z
    double LTx = (static_cast<double>(inx) / sum) * length * 2 - length;
    double LTy = -(static_cast<double>(iny) / sum * length * 2) + length;

    double RBx = (static_cast<double>(inx2 + 1) / sum) * length * 2 - length;
    double RBy = -((static_cast<double>(iny2 + 1) / sum) * length * 2) + length;

    // LT=left top, RB=right bottom
    // Returns the projected coordinates of the four corners
    std::map<std::string, std::pair<double, double>> res;
    res["LT"] = std::make_pair(LTx, LTy);
    res["RB"] = std::make_pair(RBx, RBy);
    res["LB"] = std::make_pair(LTx, RBy);
    res["RT"] = std::make_pair(RBx, LTy);

    return res;
}

// Function to get the extent based on WGS-84 coordinates
std::map<std::string, std::pair<double, double>> getExtent(double x1, double y1, double x2, double y2, int z, const std::string& source = "Google") {
    std::map<std::string, std::pair<int, int>> zb;
    
    // Calculate tile coordinates for the given WGS-84 coordinates
    int pos1x, pos1y, pos2x, pos2y;
    std::tie(pos1x, pos1y) = wgsToTile(x1, y1, z);
    std::tie(pos2x, pos2y) = wgsToTile(x2, y2, z);

    // Calculate Mercator coordinates for the four corners of the area
    std::map<std::string, std::pair<double, double>> Xframe = pixlsToMercator({
        {"LT", {pos1x, pos1y}},
        {"RT", {pos2x, pos1y}},
        {"LB", {pos1x, pos2y}},
        {"RB", {pos2x, pos2y}}
    },z);

    // Convert Mercator coordinates to WGS-84 coordinates
    for (const std::string& corner : {"LT", "LB", "RT", "RB"}) {
        std::pair<double, double>& coord = Xframe[corner];
        std::tie(coord.first, coord.second) = mercatorToWgs(coord.first, coord.second);
    }

    // Depending on the data source, apply GCJ-02 to WGS-84 conversion
    if (source == "Google") {
        // No conversion needed
    } else if (source == "Google China") {
        for (const std::string& corner : {"LT", "LB", "RT", "RB"}) {
            std::pair<double, double>& coord = Xframe[corner];
            std::tie(coord.first, coord.second) = gcj_to_wgs(coord.first, coord.second);
        }
    } else {
        throw std::invalid_argument("Invalid argument: source.");
    }

    return Xframe;
}

void get_url(std::string& url, const std::string& source, int x, int y, int z, const std::string& style) {
    if (source == "Google China") {
        url = "http://mt2.google.cn/vt/lyrs=" + style + "&hl=zh-CN&gl=CN&src=app&x=" + std::to_string(x) + "&y=" + std::to_string(y) + "&z=" + std::to_string(z);
    } else if (source == "Google") {
        url = "http://mts0.googleapis.com/vt?lyrs=" + style + "&x=" + std::to_string(x) + "&y=" + std::to_string(y) + "&z=" + std::to_string(z);
    } else {
        throw std::invalid_argument("Unknown Map Source!");
    }
}

std::vector<std::string> get_urls(double x1, double y1, double x2, double y2, int z, const std::string& source, const std::string& style) {
    std::vector<std::string> urls;
    int pos1x, pos1y, pos2x, pos2y;
    pos1x = pos1y = pos2x = pos2y = 0;

    // Get the tile coordinates
    if (source == "Google China") {
        pos1x = std::floor(x1 * std::pow(2, z));
        pos1y = std::floor(y1 * std::pow(2, z));
        pos2x = std::floor(x2 * std::pow(2, z));
        pos2y = std::floor(y2 * std::pow(2, z));
    } else if (source == "Google") {
        pos1x = std::floor((x1 + 180) * (std::pow(2, z)) / 360);
        pos1y = std::floor((1 - std::log(std::tan(y1 * M_PI / 180) + 1 / std::cos(y1 * M_PI / 180)) / M_PI) * (std::pow(2, z)) / 2);
        pos2x = std::floor((x2 + 180) * (std::pow(2, z)) / 360);
        pos2y = std::floor((1 - std::log(std::tan(y2 * M_PI / 180) + 1 / std::cos(y2 * M_PI / 180)) / M_PI) * (std::pow(2, z)) / 2);
    }

    int lenx = pos2x - pos1x + 1;
    int leny = pos2y - pos1y + 1;

    for (int j = pos1y; j < pos1y + leny; ++j) {
        for (int i = pos1x; i < pos1x + lenx; ++i) {
            std::string url;
            get_url(url, source, i, j, z, style);
            urls.push_back(url);
        }
    }
    return urls;
}

void merge_tiles(const std::vector<std::vector<char>>& datas, int x1, int y1, int x2, int y2, int z,const std::string& source, const std::string& outputFilename) {
    int pos1x, pos1y, pos2x, pos2y;
    pos1x = pos1y = pos2x = pos2y = 0;

    // Get the tile coordinates
    pos1x = std::floor(x1 * std::pow(2, z));
    pos1y = std::floor(y1 * std::pow(2, z));
    pos2x = std::floor(x2 * std::pow(2, z));
    pos2y = std::floor(y2 * std::pow(2, z));

    int lenx = pos2x - pos1x + 1;
    int leny = pos2y - pos1y + 1;
    int width = lenx * 256;
    int height = leny * 256;

    std::vector<unsigned char> mergedImage(width * height * 4, 0);

    for (size_t i = 0; i < datas.size(); ++i) {
        const char* data = datas[i].data();
        int y = i / lenx;
        int x = i % lenx;
        int xOffset = x * 256;
        int yOffset = y * 256;

        for (int py = 0; py < 256; ++py) {
            for (int px = 0; px < 256; ++px) {
                int destOffset = ((yOffset + py) * width + (xOffset + px)) * 4;
                int srcOffset = (py * 256 + px) * 4;
                mergedImage[destOffset + 0] = static_cast<unsigned char>(data[srcOffset + 2]);
                mergedImage[destOffset + 1] = static_cast<unsigned char>(data[srcOffset + 1]);
                mergedImage[destOffset + 2] = static_cast<unsigned char>(data[srcOffset + 0]);
                mergedImage[destOffset + 3] = static_cast<unsigned char>(data[srcOffset + 3]);
            }
        }
    }

    GDALAllRegister();

    GDALDataset* poDstDS;
    GDALDriver* poDriver;
    char** papszOptions = nullptr;
    poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    if (poDriver == nullptr) {
        throw std::runtime_error("GDAL driver not available.");
    }
    poDstDS = poDriver->Create(outputFilename.c_str(), width, height, 3, GDT_Byte, papszOptions);

    if (poDstDS == nullptr) {
        throw std::runtime_error("GDAL dataset creation failed.");
    }

    // 获取合并地图的四个角的空间信息并用于输出
    std::map<std::string, std::pair<double, double>> extent = getExtent(x1, y1, x2, y2, z, source);

    // 提取地理转换参数
    double geotransform[6];
    geotransform[0] = extent["LT"].first;
    geotransform[1] = (extent["RB"].first - extent["LT"].first) / width;
    geotransform[2] = 0;
    geotransform[3] = extent["LT"].second;
    geotransform[4] = 0;
    geotransform[5] = (extent["RB"].second - extent["LT"].second) / height;

    poDstDS->SetGeoTransform(geotransform);

    OGRSpatialReference oSRS;
    oSRS.importFromEPSG(4326);
    char* pszSRS_WKT = nullptr;
    oSRS.exportToWkt(&pszSRS_WKT);
    poDstDS->SetProjection(pszSRS_WKT);

    GDALRasterBand* poBand;
    for (int i = 0; i < 3; ++i) {
        poBand = poDstDS->GetRasterBand(i + 1);
        poBand->RasterIO(GF_Write, 0, 0, width, height, &mergedImage[i], width, height, GDT_Byte, 4, width * 4);
        poBand->FlushCache();
    }

    GDALClose(poDstDS);

    CPLFree(pszSRS_WKT);
}

int main() {
    if (curl_global_init(CURL_GLOBAL_DEFAULT) != CURLE_OK) {
        printf("curl_global_init failed");
        return 1;
    }

    double left = 114.43;
    double top = 30.46;
    double right = 114.46;
    double bottom = 30.43;
    int zoom = 17;
    std::string outputFilename = "output.tif";
    std::string style = "s";
    std::string source = "Google";

    // Set up GDAL environment variables
    
    std::vector<std::string> urls = get_urls(left, top, right, bottom, zoom, source, style);
    std::vector<std::vector<char>> datas(urls.size());

    int numThreads = std::thread::hardware_concurrency();
    std::vector<std::thread> threads(numThreads);
    for (int i = 0; i < numThreads; ++i) {
        threads[i] = std::thread(Downloader(i, numThreads, std::ref(urls), std::ref(datas)));
    }
    for (int i = 0; i < numThreads; ++i) {
        threads[i].join();
    }

    merge_tiles(datas, left, top, right, bottom, zoom,source, outputFilename);

    curl_global_cleanup();
    return 0;
}
