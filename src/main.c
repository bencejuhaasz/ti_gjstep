#include <tice.h>
#include <graphx.h>
#include <keypadc.h>
#include <ti/real.h>
#include <fileioc.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include <stdint.h>

/* =================== Config =================== */
#define EPS 1e-10
#define MAX_R 3               /* supports 2x3 or 3x4 augmented */
#define MAX_C (MAX_R + 1)

#define MAX_LINES 280         /* for step logs */
#define LINE_CHARS 56

/* =================== Globals =================== */
static char LOGBUF[MAX_LINES][LINE_CHARS];
static int  log_count = 0;

/* =================== Logging =================== */
static void log_line(const char *fmt, ...) {
    if (log_count >= MAX_LINES) return;
    va_list ap; va_start(ap, fmt);
    vsnprintf(LOGBUF[log_count], LINE_CHARS, fmt, ap);
    va_end(ap);
    log_count++;
}

/* =================== Fractions (smart output) =================== */
/* print "nice" fractions when possible; else compact decimal */
static void format_frac(double x, char out[LINE_CHARS]) {
    if (!isfinite(x)) { snprintf(out, LINE_CHARS, "%s", isnan(x) ? "NaN" : "inf"); return; }

    /* squash tiny noise to zero */
    if (fabs(x) < 1e-14) { snprintf(out, LINE_CHARS, "0"); return; }

    /* close to integer? */
    double rintx = round(x);
    if (fabs(x - rintx) < 1e-12) { snprintf(out, LINE_CHARS, "%.0f", rintx); return; }

    /* continued-fraction approximation with cap on denominator */
    const int64_t MAX_DEN = 1000;  /* <= 1000 keeps things readable on-screen */
    double ax = fabs(x), v = ax;
    int64_t p0=0, q0=1, p1=1, q1=0;

    for (int it=0; it<32; ++it) {
        double a = floor(v);
        int64_t p = (int64_t)a * p1 + p0;
        int64_t q = (int64_t)a * q1 + q0;

        if (q > MAX_DEN) break;

        double approx = (double)p / (double)q;
        if (fabs(approx - ax) < 5e-8) { p0=p; q0=q; p1=0; q1=0; break; }

        p0 = p1; q0 = q1; p1 = p; q1 = q;
        double r = v - a;
        if (r < 1e-15) { p0=p; q0=q; p1=0; q1=0; break; }
        v = 1.0 / r;
    }

    int have = (p1||q1);
    int64_t num = have ? p1 : p0;
    int64_t den = have ? q1 : q0;

    if (den > MAX_DEN || den == 0) {
        /* fall back to compact decimal */
        snprintf(out, LINE_CHARS, "%.6g", x);
        return;
    }

    /* put sign on numerator only */
    if (x < 0) num = -num;
    if (den < 0) { den = -den; num = -num; }

    if (den == 1) snprintf(out, LINE_CHARS, "%lld", (long long)num);
    else          snprintf(out, LINE_CHARS, "%lld/%lld", (long long)num, (long long)den);
}

static void small_val(double v, char out[16]) {
    char tmp[LINE_CHARS];
    /* zero-out ultratiny */
    if (fabs(v) < 1e-12) v = 0.0;
    format_frac(v, tmp);
    snprintf(out, 16, "%s", tmp);
}

/* =================== Homescreen input helpers =================== */
/* Parse decimal or fraction "a/b" (signs allowed). */
static bool parse_number(const char *s, double *out) {
    if (!s) return false;
    while (*s==' ' || *s=='\t') s++;
    if (!*s) { *out = 0.0; return true; }

    const char *slash = strchr(s, '/');
    if (slash) {
        char numbuf[24]={0}, denbuf[24]={0};
        size_t nlen = (size_t)(slash - s);
        if (nlen >= sizeof(numbuf)) return false;
        memcpy(numbuf, s, nlen);
        snprintf(denbuf, sizeof(denbuf), "%s", slash+1);

        double num = strtod(numbuf, NULL);
        double den = strtod(denbuf, NULL);
        if (fabs(den) < 1e-18) return false;
        *out = num / den;
        return true;
    } else {
        *out = strtod(s, NULL);
        return true;
    }
}

static double prompt_number_hs(const char *prompt) {
    char buf[32];
    while (1) {
        os_ClrHome();
        os_PutStrFull(prompt);
        memset(buf, 0, sizeof(buf));
        os_GetStringInput(NULL, buf, (uint8_t)(sizeof(buf)-1));
        double v;
        if (parse_number(buf, &v)) return v;
        os_ClrHome(); os_PutStrFull("Invalid number. Any key...");
        while (!os_GetCSC());
    }
}

static int prompt_int_hs(const char *prompt) {
    char buf[12];
    while (1) {
        os_ClrHome();
        os_PutStrFull(prompt);
        memset(buf, 0, sizeof(buf));
        os_GetStringInput(NULL, buf, (uint8_t)(sizeof(buf)-1));
        char *end = NULL;
        long v = strtol(buf, &end, 10);
        if (end != buf) return (int)v;
        os_ClrHome(); os_PutStrFull("Invalid integer. Any key...");
        while (!os_GetCSC());
    }
}

/* =================== Pretty Matrix Logger =================== */
static void log_matrix(double A[MAX_R][MAX_C], int rows, int cols) {
    log_line("Matrix [A | b]:");
    for (int i = 0; i < rows; ++i) {
        char row[LINE_CHARS]; int pos = 0;
        pos += snprintf(row+pos, LINE_CHARS-pos, "  [");
        for (int j = 0; j < cols-1; ++j) {
            char s[16]; small_val(A[i][j], s);
            pos += snprintf(row+pos, LINE_CHARS-pos, " %s", s);
        }
        char sb[16]; small_val(A[i][cols-1], sb);
        pos += snprintf(row+pos, LINE_CHARS-pos, " | %s ]", sb);
        row[LINE_CHARS-1] = 0;
        log_line("%s", row);
    }
    log_line("");
}

/* =================== Sequential input =================== */
static void sequential_input(double A[MAX_R][MAX_C], int *rows, int *cols) {
    int r = prompt_int_hs("Rows? (2 or 3): ");
    int c = prompt_int_hs("Cols? (3 or 4): ");

    if (!((r==2 && c==3) || (r==3 && c==4))) {
        os_ClrHome();
        os_PutStrFull("Only 2x3 or 3x4 allowed. Any key...");
        while (!os_GetCSC());
        r = 2; c = 3;
    }

    for (int i=0;i<r;++i) {
        for (int j=0;j<c;++j) {
            char prompt[48];
            if (j == c-1) snprintf(prompt, sizeof(prompt), "Enter b[%d]: ", i+1);
            else          snprintf(prompt, sizeof(prompt), "Enter A[%d,%d]: ", i+1, j+1);
            A[i][j] = prompt_number_hs(prompt);
        }
    }

    *rows = r; *cols = c;
}

/* =================== Gauss–Jordan (Verbose) =================== */
static void gauss_jordan_verbose(double A[MAX_R][MAX_C], int rows, int cols) {
    int iter = 1;
    log_line("Initial matrix:");
    log_matrix(A, rows, cols);

    int n = rows; /* left block is n x n */
    for (int col = 0; col < n; ++col) {
        /* pivot search */
        int pivot = col;
        double best = fabs(A[pivot][col]);
        for (int r = col + 1; r < n; ++r) {
            double v = fabs(A[r][col]);
            if (v > best) { best = v; pivot = r; }
        }
        if (best < EPS) {
            log_line("Iter %d: ~0 pivot in column %d. Singular/underdetermined.", iter++, col+1);
            log_matrix(A, rows, cols);
            return;
        }

        /* swap */
        if (pivot != col) {
            log_line("Iter %d: Swap R%d <-> R%d", iter++, col+1, pivot+1);
            for (int j=0;j<cols;++j) { double t=A[pivot][j]; A[pivot][j]=A[col][j]; A[col][j]=t; }
            log_matrix(A, rows, cols);
        }

        /* scale pivot row */
        {
            double p = A[col][col];
            if (fabs(p) < EPS) { log_line("Iter %d: pivot vanished; abort.", iter++); return; }
            double inv = 1.0 / p;

            for (int j=col;j<cols;++j) A[col][j] *= inv;

            char invs[16]; small_val(inv, invs);
            log_line("Iter %d: Scale R%d by %s (pivot->1)", iter++, col+1, invs);
            log_matrix(A, rows, cols);
        }

        /* eliminate other rows */
        for (int r=0; r<n; ++r) {
            if (r==col) continue;
            double factor = A[r][col];
            if (fabs(factor) < EPS) continue;

            for (int j=col;j<cols;++j) A[r][j] -= factor * A[col][j];

            char fs[16]; small_val(factor, fs);
            log_line("Iter %d: R%d <- R%d - (%s) * R%d", iter++, r+1, r+1, fs, col+1);
            log_matrix(A, rows, cols);
        }
    }

    log_line("Finished Gauss-Jordan. Expect [I | x].");
    log_matrix(A, rows, cols);
    log_line("Solution x:");
    for (int i=0;i<rows;++i){ char s[16]; small_val(A[i][cols-1], s); log_line("  x[%d] = %s", i, s); }
}

/* =================== GraphX scroll viewer =================== */
static void show_log_viewer(void) {
    const int margin = 4;
    const int line_h = 8;
    const int lines_on_screen = (LCD_HEIGHT - 2*margin) / line_h;
    int top = 0;

    gfx_Begin();
    gfx_SetDrawBuffer();

    for (;;) {
        kb_Scan();
        if (kb_Data[6] & kb_Clear) break;
        if (kb_Data[7] & kb_Up)    { if (top > 0) top--; delay(16); }
        if (kb_Data[7] & kb_Down)  { if (top + lines_on_screen < log_count) top++; delay(16); }
        if (kb_Data[7] & kb_Left)  { top -= lines_on_screen; if (top < 0) top = 0; delay(60); }
        if (kb_Data[7] & kb_Right) { top += lines_on_screen; if (top > log_count - lines_on_screen) top = log_count - lines_on_screen; if (top < 0) top = 0; delay(60); }

        gfx_FillScreen(255);
        gfx_SetTextFGColor(0);
        gfx_SetTextBGColor(255);
        gfx_SetTextScale(1,1);

        gfx_PrintStringXY("Gauss-Jordan Steps (UP/DOWN, CLEAR exit)", margin, margin);

        int y = margin + line_h + 2;
        int shown = 0;
        for (int i = top; i < log_count && shown < lines_on_screen-2; ++i, ++shown) {
            gfx_PrintStringXY(LOGBUF[i], margin, y);
            y += line_h;
        }

        char footer[60];
        snprintf(footer, sizeof(footer), "Lines %d-%d / %d", top+1, top+shown, log_count);
        gfx_PrintStringXY(footer, margin, LCD_HEIGHT - margin - line_h);

        gfx_SwapDraw();
    }

    gfx_End();
}

/* =================== main =================== */
int main(void) {
    double A[MAX_R][MAX_C] = {{0}};
    int rows=0, cols=0;

    sequential_input(A, &rows, &cols);

    log_count = 0;
    gauss_jordan_verbose(A, rows, cols);
    show_log_viewer();
    return 0;
}
