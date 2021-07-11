#include <bits/stdc++.h>
using namespace std;

class Graphics {
public:
	double screenW;
	double screenH;
	ostringstream data;

	double sr;
	double sg;
	double sb;
	double sa;
	double fr;
	double fg;
	double fb;
	double fa;

	Graphics() : screenW(1), screenH(1), sr(0), sg(0), sb(0), sa(1), fr(1), fg(1), fb(1), fa(1) {
	}

	void screen(int width, int height) {
		screenW = width;
		screenH = height;
	}

	void clear() {
		data.str("");
		data.clear(stringstream::goodbit);
	}

	void stroke(double r, double g, double b) {
		stroke(r, g, b, 1);
	}

	void stroke(double r, double g, double b, double a) {
		sr = r;
		sg = g;
		sb = b;
		sa = a;
	}

	void noStroke() {
		stroke(0, 0, 0, 0);
	}

	void fill(double r, double g, double b) {
		fill(r, g, b, 1);
	}

	void fill(double r, double g, double b, double a) {
		fr = r;
		fg = g;
		fb = b;
		fa = a;
	}

	void noFill() {
		fill(0, 0, 0, 0);
	}

	void line(double x1, double y1, double x2, double y2) {
		data << "<line x1=\"" << x1 << "\" y1=\"" << y1 << "\" x2=\"" << x2 << "\" y2=\"" << y2 << "\" " << stroke() << "/>\n";
	}

	void rect(double x, double y, double w, double h) {
		data << "<rect x=\"" << x << "\" y=\"" << y << "\" width=\"" << w << "\" height=\"" << h << "\" " << stroke() << " " + fill() << "/>\n";
	}

    void poly(vector<double> x, vector<double> y) {
        assert(x.size() == y.size());
        data << "<polygon points=\"";
        for(size_t i=0; i<x.size(); i++) {
            data << x[i] << "," << y[i];
            if(i + 1 < x.size()) data << " ";
        }
        data << "\" " << stroke() << " " << fill() << "/>\n";
    }

    void poly(vector< complex<double> > pts) {
        vector<double> xs, ys;
        for(auto pt : pts) {
            double x = pt.real(), y = pt.imag();
            xs.emplace_back(x);
            ys.emplace_back(y);
        }
        poly(xs, ys);
    }

	void text(string str, double x, double y, double size = 16) {
		data << "<text text-anchor=\"middle\" x=\"" << x << "\" y=\"" << y << "\" font-size=\"" << size << "\" " << fill() << " >" << str << "</text>\n";
	}

	string dump(string id = "", string style = "", int widthPx = -1, int heightPx = -1) const {
		ostringstream res;
		res << "<svg ";
		if (id != "") res << "id=\"" + id + "\" ";
		if (style != "") res << "style=\"" + style + "\" ";
		if (widthPx != -1) res << "width=\"" << widthPx << "\" ";
		if (heightPx != -1) res << "height=\"" << heightPx << "\" ";
		res << "viewBox=\"-1 -1 " << screenW + 2 << " " << screenH + 2 << "\" xmlns=\"http://www.w3.org/2000/svg\">\n" << data.str() << "</svg>";
		return res.str();
	}

private:
	string stroke() const {
		return "stroke=\"" + rgb(sr, sg, sb) + "\" stroke-opacity=\"" + s(sa) + "\"";
	}

	string fill() const {
		return "fill=\"" + rgb(fr, fg, fb) + "\" fill-opacity=\"" + s(fa) + "\"";
	}

	string rgb(double r, double g, double b) const {
		return "rgb(" + s(lround(r * 255)) + "," + s(lround(g * 255)) + "," + s(lround(b * 255)) + ")";
	}

	string s(double x) const {
		return to_string(x);
	}
};

class Movie {
public:
	vector<string> svgs;

	Movie() {
	}

	void clear() {
		svgs.clear();
	}

	void addFrame(Graphics& g) {
		svgs.push_back(g.dump("f" + to_string(svgs.size()), "display:none;pointer-events:none;user-select:none;"));
	}

	string dumpHtml(int fps) {
		ostringstream res;
		res << "<html><body><div id=\"text\">loading...</div>" << endl;
		for (string& svg : svgs) {
			res << svg << endl;
		}

		res << "<script>\nlet numFrames = " << svgs.size() << ", fps = " << fps << ";";
		res << R"(
	let text = document.getElementById("text");
	let frames = [];
	for (let i = 0; i < numFrames; i++) {
		let f = document.getElementById("f" + i);
		frames.push(f);
		f.style.display = "none";
	}
	let currentFrame = 0;
	let playing = true;
	setInterval(() => {
		if (!playing) return;
		text.innerText = (currentFrame + 1) + " / " + numFrames;
		frames[(currentFrame - 1 + numFrames) % numFrames].style.display = "none";
		frames[currentFrame].style.display = null;
		currentFrame = (currentFrame + 1) % numFrames;
		if (currentFrame == 0) playing = false;
	}, 1000 / fps);
	window.onmousedown = e => { if (e.button == 0) playing = true; };
;)";
		res << "</script>" << endl;
		res << "</body></html>" << endl;
		return res.str();
	}
private:
};
