#if defined(__ARM_NEON)
	#include <arm_neon.h>
#endif
#include <algorithm>
#include <memory>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdint.h>

enum
{
	// get
	kMethodGet1,
	kMethodGet4,
	kMethodGet8,
	kMethodGet16,

	// set
	kMethodMemset,
	kMethodSet1,
	kMethodSet4,
	kMethodSet8,
	kMethodSet16,

	// reverse
	kMethodStdReverse,
	kMethodReverse1,
	kMethodReverse4,
	kMethodReverse8,
	kMethodReverse16,

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
	uint8_t exec(volatile uint8_t* p, uint8_t v, size_t size) const override { return 0; }
};

template<> struct Method<kMethodGet1> : public IMethod
{
	const char* name() const override { return "get1"; }

	__attribute__((noinline))
	uint8_t exec(volatile uint8_t* p, uint8_t v, size_t size) const override
	{
		for ( size_t i = 0; i < size; i++ ) {
			v ^= p[i];
		}
		return v;
	}
};

template<> struct Method<kMethodGet4> : public IMethod
{
	const char* name() const override { return "get4"; }

	__attribute__((noinline))
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

template<> struct Method<kMethodGet8> : public IMethod
{
	const char* name() const override { return "get8"; }

	__attribute__((noinline))
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

#if defined(__ARM_NEON)
template<> struct Method<kMethodGet16> : public IMethod
{
	const char* name() const override { return "get16"; }

	__attribute__((noinline))
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

template<> struct Method<kMethodMemset> : public IMethod
{
	const char* name() const override { return "memset"; }

	__attribute__((noinline))
	uint8_t exec(volatile uint8_t* p, uint8_t v, size_t size) const override
	{
		std::memset(const_cast<uint8_t*>(p), v, size);
		return v;
	}
};

template<> struct Method<kMethodSet1> : public IMethod
{
	const char* name() const override { return "set1"; }

	__attribute__((noinline))
	uint8_t exec(volatile uint8_t* p, uint8_t v, size_t size) const override
	{
		for ( size_t i = 0; i < size; i++ ) {
			p[i] = v;
		}
		return v;
	}
};

template<> struct Method<kMethodSet4> : public IMethod
{
	const char* name() const override { return "set4"; }

	__attribute__((noinline))
	uint8_t exec(volatile uint8_t* p, uint8_t v, size_t size) const override
	{
		uint32_t w = v;
		for ( size_t i = 0; i < size; i += sizeof(w) ) {
			*reinterpret_cast<volatile uint32_t*>(&p[i]) = w;
		}
		return v;
	}
};

template<> struct Method<kMethodSet8> : public IMethod
{
	const char* name() const override { return "set8"; }

	__attribute__((noinline))
	uint8_t exec(volatile uint8_t* p, uint8_t v, size_t size) const override
	{
		uint64_t w = v;
		for ( size_t i = 0; i < size; i += sizeof(w) ) {
			*reinterpret_cast<volatile uint64_t*>(&p[i]) = w;
		}
		return v;
	}
};

#if defined(__ARM_NEON)
template<> struct Method<kMethodSet16> : public IMethod
{
	const char* name() const override { return "set16"; }

	__attribute__((noinline))
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

template<> struct Method<kMethodStdReverse> : public IMethod
{
	const char* name() const override { return "std::reverse"; }

	__attribute__((noinline))
	uint8_t exec(volatile uint8_t* p, uint8_t v, size_t size) const override
	{
		std::reverse(p, &p[size - 1]);
		return v;
	}
};

template<> struct Method<kMethodReverse1> : public IMethod
{
	const char* name() const override { return "reverse1"; }

	__attribute__((noinline))
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

template<> struct Method<kMethodReverse4> : public IMethod
{
	const char* name() const override { return "reverse4"; }

	__attribute__((noinline))
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

template<> struct Method<kMethodReverse8> : public IMethod
{
	const char* name() const override { return "reverse8"; }

	__attribute__((noinline))
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

#if defined(__ARM_NEON)
template<> struct Method<kMethodReverse16> : public IMethod
{
	const char* name() const override { return "reverse16"; }

	__attribute__((noinline))
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

static IMethod *g_Methods[] = {
	new Method<kMethodGet1>,
	new Method<kMethodGet4>,
	new Method<kMethodGet8>,
	new Method<kMethodGet16>,
	new Method<kMethodMemset>,
	new Method<kMethodSet1>,
	new Method<kMethodSet4>,
	new Method<kMethodSet8>,
	new Method<kMethodSet16>,
	new Method<kMethodStdReverse>,
	new Method<kMethodReverse1>,
	new Method<kMethodReverse4>,
	new Method<kMethodReverse8>,
	new Method<kMethodReverse16>,
};

static const size_t kAlignment = 4096;

static int Randomize(uint8_t* p, size_t size)
{
	uint8_t random[kAlignment];
	FILE* urandomFile = fopen("/dev/urandom", "rb");
	size_t numRead = fread(random, 1, sizeof(random), urandomFile);
	fclose(urandomFile);
	if ( sizeof(random) != numRead ) {
		std::printf("numRead: %lu\n", numRead);
		return 1;
	}

	for ( size_t i = 0, j = 0; i < size; i += kAlignment/2, j++ ) {
		memcpy(&p[i], &random[j & (kAlignment-1)], kAlignment/2);
	}

	return 0;
}

static void Usage(const char* argv0)
{
	std::printf("usage\n");
	std::printf("\t%s cycles size method\n\n", argv0);
	std::printf("cycles\n");
	std::printf("\tNumber of times to scan memory.\n");
	std::printf("size\n");
	std::printf("\tSize of memory allocation. (align up to %lu)\n", kAlignment);
	std::printf("method\n");
	std::printf("\tThere are kind of 3 type functions that memory read, memory write and memory read-write.\n");
	std::printf("\tAnd, there are type of element widths. For example, get1 accesses 1-byte at a time.\n");
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

	if ( int err = Randomize(p, size) ) {
		return err;
	}

	uint8_t v = uint8_t(numCycles);
	for ( int i = 0; i < numCycles; i++ ) {
		v = method->exec(p, v, size);
	}
	std::printf("result (dummy data): %d\n", v);

	return 0;
}
