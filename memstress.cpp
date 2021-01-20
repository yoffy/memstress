#if defined(__SSE2__) || defined(__AVX2__)
	#include <x86intrin.h>
#endif
#if defined(__ARM_NEON)
	#include <arm_neon.h>
#endif
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <random>
#include <stdint.h>

enum
{
	// read
	kMethodRead1,
	kMethodRead4,
	kMethodRead8,
	kMethodRead16,
	kMethodRead32,

	// write
	kMethodStdMemset,
	kMethodWrite1,
	kMethodWrite4,
	kMethodWrite8,
	kMethodWrite16,
	kMethodWrite32,

	// readwrite
	kMethodStdReverse,
	kMethodReadWrite1,
	kMethodReadWrite4,
	kMethodReadWrite8,
	kMethodReadWrite16,
	kMethodReadWrite32,

	kMethodEnd
};

struct IMethod
{
	virtual const char* name() const = 0;
	virtual uint8_t exec(volatile uint8_t* p, uint8_t v, size_t size) const = 0;
};

template<int N>
struct Method : public IMethod
{
	const char* name() const override { return nullptr; }
	uint8_t exec(volatile uint8_t*, uint8_t, size_t) const override { return 0; }
};

template<> struct Method<kMethodRead1> : public IMethod
{
	const char* name() const override { return "read1"; }

	uint8_t exec(volatile uint8_t* p, uint8_t v, size_t size) const override
	{
		for ( size_t i = 0; i < size; i++ ) {
			v ^= p[i];
		}
		return v;
	}
};

template<> struct Method<kMethodRead4> : public IMethod
{
	const char* name() const override { return "read4"; }

	uint8_t exec(volatile uint8_t* p, uint8_t v, size_t size) const override
	{
		uint32_t w = 0;
		for ( size_t i = 0; i < size; i += sizeof(w) ) {
			w ^= *reinterpret_cast<const volatile uint32_t*>(&p[i]);
		}

		v = v
			^ ( w        & 0xFF)
			^ ((w >>  8) & 0xFF)
			^ ((w >> 16) & 0xFF)
			^ ((w >> 24) & 0xFF);
		return v;
	}
};

template<> struct Method<kMethodRead8> : public IMethod
{
	const char* name() const override { return "read8"; }

	uint8_t exec(volatile uint8_t* p, uint8_t v, size_t size) const override
	{
		uint64_t w = 0;
		for ( size_t i = 0; i < size; i += sizeof(w) ) {
			w ^= *reinterpret_cast<const volatile uint64_t*>(&p[i]);
		}

		v = v
			^ ( w        & 0xFF)
			^ ((w >>  8) & 0xFF)
			^ ((w >> 16) & 0xFF)
			^ ((w >> 24) & 0xFF)
			^ ((w >> 32) & 0xFF)
			^ ((w >> 40) & 0xFF)
			^ ((w >> 48) & 0xFF)
			^ ((w >> 56) & 0xFF);
		return v;
	}
};

#if defined(__SSE2__)
template<> struct Method<kMethodRead16> : public IMethod
{
	const char* name() const override { return "read16"; }

	uint8_t exec(volatile uint8_t* p, uint8_t v, size_t size) const override
	{
		__m128i w = _mm_setzero_si128();
		for ( size_t i = 0; i < size; i += sizeof(w) ) {
			w ^= _mm_loadu_si128((const __m128i*)(&p[i]));
		}

		uint64_t x = _mm_cvtsi128_si64(w) ^ _mm_cvtsi128_si64(_mm_unpackhi_epi64(w, w));
		v = v
			^ ( x        & 0xFF)
			^ ((x >>  8) & 0xFF)
			^ ((x >> 16) & 0xFF)
			^ ((x >> 24) & 0xFF)
			^ ((x >> 32) & 0xFF)
			^ ((x >> 40) & 0xFF)
			^ ((x >> 48) & 0xFF)
			^ ((x >> 56) & 0xFF);
		return v;
	}
};
#endif

#if defined(__ARM_NEON)
template<> struct Method<kMethodRead16> : public IMethod
{
	const char* name() const override { return "read16"; }

	uint8_t exec(volatile uint8_t* p, uint8_t v, size_t size) const override
	{
		uint8x16_t w = vdupq_n_u8(0);
		for ( size_t i = 0; i < size; i += sizeof(w) ) {
			w ^= vld1q_u8(const_cast<const uint8_t*>(&p[i]));
		}

		uint64_t x = vgetq_lane_u64(vreinterpretq_u64_u8(w), 0) ^ vgetq_lane_u64(vreinterpretq_u64_u8(w), 1);
		v = v
			^ ( x        & 0xFF)
			^ ((x >>  8) & 0xFF)
			^ ((x >> 16) & 0xFF)
			^ ((x >> 24) & 0xFF)
			^ ((x >> 32) & 0xFF)
			^ ((x >> 40) & 0xFF)
			^ ((x >> 48) & 0xFF)
			^ ((x >> 56) & 0xFF);
		return v;
	}
};
#endif

#if defined(__AVX2__)
template<> struct Method<kMethodRead32> : public IMethod
{
	const char* name() const override { return "read32"; }

	uint8_t exec(volatile uint8_t* p, uint8_t v, size_t size) const override
	{
		__m256i w = _mm256_setzero_si256();
		for ( size_t i = 0; i < size; i += sizeof(w) ) {
			w ^= _mm256_loadu_si256((const __m256i*)(&p[i]));
		}

		__m128i w128 = _mm256_castsi256_si128(w) ^ _mm256_extracti128_si256(w, 1);
		uint64_t x = _mm_cvtsi128_si64(w128) ^ _mm_cvtsi128_si64(_mm_unpackhi_epi64(w128, w128));
		v = v
			^ ( x        & 0xFF)
			^ ((x >>  8) & 0xFF)
			^ ((x >> 16) & 0xFF)
			^ ((x >> 24) & 0xFF)
			^ ((x >> 32) & 0xFF)
			^ ((x >> 40) & 0xFF)
			^ ((x >> 48) & 0xFF)
			^ ((x >> 56) & 0xFF);
		return v;
	}
};
#endif

template<> struct Method<kMethodStdMemset> : public IMethod
{
	const char* name() const override { return "std::memset"; }

	uint8_t exec(volatile uint8_t* p, uint8_t v, size_t size) const override
	{
		std::memset(const_cast<uint8_t*>(p), v, size);
		return v;
	}
};

template<> struct Method<kMethodWrite1> : public IMethod
{
	const char* name() const override { return "write1"; }

	uint8_t exec(volatile uint8_t* p, uint8_t v, size_t size) const override
	{
		for ( size_t i = 0; i < size; i++ ) {
			p[i] = v;
		}
		return v;
	}
};

template<> struct Method<kMethodWrite4> : public IMethod
{
	const char* name() const override { return "write4"; }

	uint8_t exec(volatile uint8_t* p, uint8_t v, size_t size) const override
	{
		uint32_t w = v;
		for ( size_t i = 0; i < size; i += sizeof(w) ) {
			*reinterpret_cast<volatile uint32_t*>(&p[i]) = w;
		}
		return v;
	}
};

template<> struct Method<kMethodWrite8> : public IMethod
{
	const char* name() const override { return "write8"; }

	uint8_t exec(volatile uint8_t* p, uint8_t v, size_t size) const override
	{
		uint64_t w = v;
		for ( size_t i = 0; i < size; i += sizeof(w) ) {
			*reinterpret_cast<volatile uint64_t*>(&p[i]) = w;
		}
		return v;
	}
};

#if defined(__SSE2__)
template<> struct Method<kMethodWrite16> : public IMethod
{
	const char* name() const override { return "write16"; }

	uint8_t exec(volatile uint8_t* p, uint8_t v, size_t size) const override
	{
		__m128i w = _mm_set1_epi8(v);
		for ( size_t i = 0; i < size; i += sizeof(w) ) {
			_mm_storeu_si128((__m128i*)(&p[i]), w);
		}
		return v;
	}
};
#endif

#if defined(__ARM_NEON)
template<> struct Method<kMethodWrite16> : public IMethod
{
	const char* name() const override { return "write16"; }

	uint8_t exec(volatile uint8_t* p, uint8_t v, size_t size) const override
	{
		uint8x16_t w = vdupq_n_u8(v);
		for ( size_t i = 0; i < size; i += sizeof(w) ) {
			vst1q_u8(const_cast<uint8_t*>(&p[i]), w);
		}
		return v;
	}
};
#endif

#if defined(__AVX2__)
template<> struct Method<kMethodWrite32> : public IMethod
{
	const char* name() const override { return "write32"; }

	uint8_t exec(volatile uint8_t* p, uint8_t v, size_t size) const override
	{
		__m256i w = _mm256_set1_epi8(v);
		for ( size_t i = 0; i < size; i += sizeof(w) ) {
			_mm256_storeu_si256((__m256i*)(&p[i]), w);
		}
		return v;
	}
};
#endif

template<> struct Method<kMethodStdReverse> : public IMethod
{
	const char* name() const override { return "std::reverse"; }

	uint8_t exec(volatile uint8_t* p, uint8_t v, size_t size) const override
	{
		std::reverse(p, &p[size - 1]);
		return v;
	}
};

template<> struct Method<kMethodReadWrite1> : public IMethod
{
	const char* name() const override { return "readwrite1"; }

	uint8_t exec(volatile uint8_t* p, uint8_t v, size_t size) const override
	{
		for ( size_t i = 0; i < size/2; i++ ) {
			size_t j = size - 1 - i;
			uint8_t tmp1 = p[i];
			uint8_t tmp2 = p[j];
			p[i] = tmp2;
			p[j] = tmp1;
		}
		return v;
	}
};

template<> struct Method<kMethodReadWrite4> : public IMethod
{
	const char* name() const override { return "readwrite4"; }

	uint8_t exec(volatile uint8_t* p, uint8_t v, size_t size) const override
	{
		for ( size_t i = 0; i < size/2; i += sizeof(uint32_t) ) {
			size_t j = size - sizeof(uint32_t) - i;
			uint32_t tmp1 = *reinterpret_cast<volatile uint32_t*>(&p[i]);
			uint32_t tmp2 = *reinterpret_cast<volatile uint32_t*>(&p[j]);
			*reinterpret_cast<volatile uint32_t*>(&p[i]) = byterev(tmp2);
			*reinterpret_cast<volatile uint32_t*>(&p[j]) = byterev(tmp1);
		}
		return v;
	}

private:
	static uint32_t byterev(uint32_t v)
	{
		return ((v        & 0xFF) << 24)
			| (((v >> 8 ) & 0xFF) << 16)
			| (((v >> 16) & 0xFF) <<  8)
			| (((v >> 24) & 0xFF)      );
	}
};

template<> struct Method<kMethodReadWrite8> : public IMethod
{
	const char* name() const override { return "readwrite8"; }

	uint8_t exec(volatile uint8_t* p, uint8_t v, size_t size) const override
	{
		for ( size_t i = 0; i < size/2; i += sizeof(uint64_t) ) {
			size_t j = size - sizeof(uint64_t) - i;
			uint64_t tmp1 = *reinterpret_cast<volatile uint64_t*>(&p[i]);
			uint64_t tmp2 = *reinterpret_cast<volatile uint64_t*>(&p[j]);
			*reinterpret_cast<volatile uint64_t*>(&p[i]) = byterev(tmp2);
			*reinterpret_cast<volatile uint64_t*>(&p[j]) = byterev(tmp1);
		}
		return v;
	}

private:
	static uint64_t byterev(uint64_t v)
	{
		return ((v        & 0xFF) << 56)
			| (((v >> 8 ) & 0xFF) << 48)
			| (((v >> 16) & 0xFF) << 40)
			| (((v >> 24) & 0xFF) << 32)
			| (((v >> 32) & 0xFF) << 24)
			| (((v >> 40) & 0xFF) << 16)
			| (((v >> 48) & 0xFF) <<  8)
			| (((v >> 56) & 0xFF)      );
	}
};

#if defined(__SSSE3__)
template<> struct Method<kMethodReadWrite16> : public IMethod
{
	const char* name() const override { return "readwrite16"; }

	uint8_t exec(volatile uint8_t* p, uint8_t v, size_t size) const override
	{
		for ( size_t i = 0; i < size/2; i += sizeof(__m128i) ) {
			size_t j = size - sizeof(__m128i) - i;
			__m128i tmp1 = _mm_loadu_si128((const __m128i*)(&p[i]));
			__m128i tmp2 = _mm_loadu_si128((const __m128i*)(&p[j]));
			_mm_storeu_si128((__m128i*)(&p[i]), byterev(tmp2));
			_mm_storeu_si128((__m128i*)(&p[j]), byterev(tmp1));
		}
		return v;
	}

private:
	static __m128i byterev(__m128i v)
	{
		__m128i table = _mm_set_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
		return _mm_shuffle_epi8(v, table);
	}
};
#endif

#if defined(__ARM_NEON)
template<> struct Method<kMethodReadWrite16> : public IMethod
{
	const char* name() const override { return "readwrite16"; }

	uint8_t exec(volatile uint8_t* p, uint8_t v, size_t size) const override
	{
		for ( size_t i = 0; i < size/2; i += sizeof(uint8x16_t) ) {
			size_t j = size - sizeof(uint8x16_t) - i;
			uint8x16_t tmp1 = vld1q_u8(const_cast<const uint8_t*>(&p[i]));
			uint8x16_t tmp2 = vld1q_u8(const_cast<const uint8_t*>(&p[j]));
			vst1q_u8(const_cast<uint8_t*>(&p[i]), byterev(tmp2));
			vst1q_u8(const_cast<uint8_t*>(&p[j]), byterev(tmp1));
		}
		return v;
	}

private:
	static uint8x16_t byterev(uint8x16_t v)
	{
		// 32 10 -> 23 01
		uint8x16_t rev64 = vrev64q_u8(v);
		// 23 01 -> 01 23
		return vcombine_u8(vget_high_u8(rev64), vget_low_u8(rev64));
	}
};
#endif

#if defined(__AVX2__)
template<> struct Method<kMethodReadWrite32> : public IMethod
{
	const char* name() const override { return "readwrite32"; }

	uint8_t exec(volatile uint8_t* p, uint8_t v, size_t size) const override
	{
		for ( size_t i = 0; i < size/2; i += sizeof(__m256i) ) {
			size_t j = size - sizeof(__m256i) - i;
			__m256i tmp1 = _mm256_loadu_si256((const __m256i*)(&p[i]));
			__m256i tmp2 = _mm256_loadu_si256((const __m256i*)(&p[j]));
			_mm256_storeu_si256((__m256i*)(&p[i]), byterev(tmp2));
			_mm256_storeu_si256((__m256i*)(&p[j]), byterev(tmp1));
		}
		return v;
	}

private:
	static __m256i byterev(__m256i v)
	{
		__m256i table = _mm256_set_epi8(
			0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
			0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15
		);
		// 32 10 -> 23 01
		__m256i rev128 = _mm256_shuffle_epi8(v, table);
		// 23 01 -> 01 23
		__m128i lo = _mm256_castsi256_si128(rev128);
		__m128i hi = _mm256_extracti128_si256(rev128, 1);
		return _mm256_inserti128_si256(_mm256_castsi128_si256(hi), lo, 1);
	}
};
#endif

static IMethod *g_Methods[] = {
	new Method<kMethodRead1>,
	new Method<kMethodRead4>,
	new Method<kMethodRead8>,
	new Method<kMethodRead16>,
	new Method<kMethodRead32>,
	new Method<kMethodStdMemset>,
	new Method<kMethodWrite1>,
	new Method<kMethodWrite4>,
	new Method<kMethodWrite8>,
	new Method<kMethodWrite16>,
	new Method<kMethodWrite32>,
	new Method<kMethodStdReverse>,
	new Method<kMethodReadWrite1>,
	new Method<kMethodReadWrite4>,
	new Method<kMethodReadWrite8>,
	new Method<kMethodReadWrite16>,
	new Method<kMethodReadWrite32>,
};

static const size_t kAlignment = 4096;

static uint64_t Random()
{
	static uint64_t x = std::random_device()() | (1ULL << 63);
	x = x ^ (x << 7);
	return x = x ^ (x >> 9);
}

static int Initialize(uint8_t* p, size_t size)
{
	for ( size_t i = 0; i < size; i += sizeof(uint64_t) ) {
		*reinterpret_cast<uint64_t*>(&p[i]) = Random();
	}
	return 0;
}

static void Usage(const char* argv0)
{
	std::printf("usage\n");
	std::printf("\t%s cycles size method\n\n", argv0);
	std::printf("cycles\n");
	std::printf("\tNumber of times to scan memory. (0 is infinite)\n");
	std::printf("size\n");
	std::printf("\tSize of memory allocation. (align up to %lu)\n", kAlignment);
	std::printf("method\n");
	std::printf("\tThere are kind of 3 type functions that memory read, memory write and memory read-write.\n");
	std::printf("\tAnd, there are type of element widths. For example, read1 accesses 1-byte at a time.\n");
	std::printf("\tUse following methods:\n");
	std::printf("\t\t");
	int numMethods = sizeof(g_Methods) / sizeof(*g_Methods);
	for ( int i = 0; i < numMethods; i++ ) {
		if ( const char* name = g_Methods[i]->name() ) {
			std::printf("%s ", name);
		}
	}
	std::printf("\n");
}

static size_t AlignUp(size_t v, size_t align)
{
	return -(-v & -align);
}

int main(int argc, char* argv[])
{
	if ( argc != 4 ) {
		Usage(argv[0]);
		return 1;
	}

	int numCycles = std::atoi(argv[1]);
	long size = AlignUp(std::atol(argv[2]), kAlignment);
	const char* cstrMethod = argv[3];

	IMethod* method = nullptr;
	int numMethods = sizeof(g_Methods) / sizeof(*g_Methods);
	for ( int i = 0; i < numMethods; i++ ) {
		if ( const char* name = g_Methods[i]->name() ) {
			if ( std::strcmp(name, cstrMethod) == 0 ) {
				method = g_Methods[i];
			}
		}
	}
	if ( ! method ) {
		Usage(argv[0]);
		std::printf("\nInvalid method.\n");
		return 1;
	}
	std::printf("cycles: %d\n", numCycles);
	std::printf("size: %ld\n", size);
	std::printf("method: %s\n", method->name());

	uint8_t* p = reinterpret_cast<uint8_t*>(malloc(size));

	if ( int err = Initialize(p, size) ) {
		return err;
	}

	uint8_t v = uint8_t(numCycles);
	for ( int i = 0; ! numCycles || i < numCycles; i++ ) {
		v = method->exec(p, v, size);
	}
	std::printf("result (dummy data): %d\n", v);

	return 0;
}
