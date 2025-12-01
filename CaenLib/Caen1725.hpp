#pragma once

#include "CAENDigitizer.h"
#include <string>
#include <stdexcept>
#include <vector>
#include <cstdio>    // snprintf
#include <utility>   // std::exchange

class Digitizer {
    int handle_ = -1;

public:
    // ------------------------------------------------------------------
    // Construction / destruction
    // ------------------------------------------------------------------
    Digitizer() = default;

    // Delete copy, allow move
    Digitizer(const Digitizer&) = delete;
    Digitizer& operator=(const Digitizer&) = delete;

    Digitizer(Digitizer&& other) noexcept : handle_(std::exchange(other.handle_, -1)) {}
    Digitizer& operator=(Digitizer&& other) noexcept {
        if (this != &other) {
            close();
            handle_ = std::exchange(other.handle_, -1);
        }
        return *this;
    }

    ~Digitizer() noexcept { close(); }

    // ------------------------------------------------------------------
    // Connect – tries USB first, then OpticalLink
    // ------------------------------------------------------------------
    // linkNum: USB device number or optical link number
    // conetNode: daisy-chain node (for optical) or USB device number for USB (same as linkNum)
    void connect(int linkNum = 0, int conetNode = 0, uint32_t vmeBase = 0)
    {
        close();  // safety

        // Try USB first (arg must be nullptr; linkNum -> ConetNode parameter)
        CAEN_DGTZ_ErrorCode ret = CAEN_DGTZ_OpenDigitizer2(
            CAEN_DGTZ_USB,
            nullptr,
            linkNum,      // goes to ConetNode parameter for USB device number
            0,
            &handle_);

        if (ret != CAEN_DGTZ_Success) {
            // Try Optical Link (common case: arg == nullptr, conetNode is the daisy node)
            ret = CAEN_DGTZ_OpenDigitizer2(
                CAEN_DGTZ_OpticalLink,
                nullptr,
                conetNode,
                vmeBase,
                &handle_);
        }

        if (ret != CAEN_DGTZ_Success || handle_ < 0) {
            handle_ = -1;
            throw std::runtime_error("Cannot open digitizer (USB nor OpticalLink)");
        }

        // Optionally reset to known default state (safer). Remove if you must preserve board config.
        cmd(CAEN_DGTZ_Reset(handle_), "Reset after Open");
    }

    // ------------------------------------------------------------------
    // Explicit connect with exact link type
    // ------------------------------------------------------------------
    void connect(CAEN_DGTZ_ConnectionType type,
             int linkNum = 0, int conetNode = 0, uint32_t vmeBase = 0)
    {
        close();

        const void* arg = nullptr;
        int conetNodeParam = 0;
        uint32_t vmeParam = vmeBase;

        switch (type) {
            case CAEN_DGTZ_USB:
                // USB: linkNum goes into ConetNode parameter
                conetNodeParam = linkNum;
                vmeParam = 0;
                break;

            case CAEN_DGTZ_OpticalLink:
                // Optical: conetNode used for daisy-chain
                conetNodeParam = conetNode;
                break;

            case CAEN_DGTZ_USB_A4818:
                // A4818 USB: same logic as USB
                conetNodeParam = linkNum;
                break;

            default:
                // Unknown or unsupported types:
                // We pass through conetNodeParam and vmeParam untouched.
                conetNodeParam = conetNode;
                break;
        }

        CAEN_DGTZ_ErrorCode ret =
            CAEN_DGTZ_OpenDigitizer2(type, arg, conetNodeParam, vmeParam, &handle_);

        if (ret != CAEN_DGTZ_Success || handle_ < 0) {
            handle_ = -1;
            throw std::runtime_error("Failed to open digitizer with given link type");
        }

        cmd(CAEN_DGTZ_Reset(handle_), "Reset after Open");
    }

    // ------------------------------------------------------------------
    // Basic operations – thin wrappers that throw on error
    // ------------------------------------------------------------------
    void cmd(CAEN_DGTZ_ErrorCode ret, const char* where) const
    {
        if (ret != CAEN_DGTZ_Success) {
            char msg[256];
            snprintf(msg, sizeof(msg), "CAEN error %d at %s", static_cast<int>(ret), where);
            throw std::runtime_error(msg);
        }
    }

    // convenience wrapper keeping previous single-arg calls (reports "cmd")
    void cmd(CAEN_DGTZ_ErrorCode ret) const { cmd(ret, "CAEN API call"); }

    void reset()            { cmd(CAEN_DGTZ_Reset(handle_), "Reset"); }
    void start()            { cmd(CAEN_DGTZ_SWStartAcquisition(handle_), "SWStartAcquisition"); }
    void stop()             { cmd(CAEN_DGTZ_SWStopAcquisition(handle_), "SWStopAcquisition"); }
    void softwareTrigger()  { cmd(CAEN_DGTZ_SendSWtrigger(handle_), "SendSWtrigger"); }

    void setRecordLength(uint32_t samples)        { cmd(CAEN_DGTZ_SetRecordLength(handle_, samples), "SetRecordLength"); }
    void setChannelMask(uint32_t mask)            { cmd(CAEN_DGTZ_SetChannelEnableMask(handle_, mask), "SetChannelEnableMask"); }
    void setDCOffset(uint32_t ch, uint32_t value) { cmd(CAEN_DGTZ_SetChannelDCOffset(handle_, ch, value), "SetChannelDCOffset"); }
    void setPostTrigger(uint32_t percent)         { cmd(CAEN_DGTZ_SetPostTriggerSize(handle_, percent), "SetPostTriggerSize"); }

    // ------------------------------------------------------------------
    // Readout buffer – managed automatically per board
    // ------------------------------------------------------------------
    class Buffer {
        char* ptr = nullptr;
        uint32_t size = 0;
        int handle = -1;
        friend class Digitizer;
    public:
        Buffer() = default;
        explicit Buffer(int h) : handle(h) {
            // CAEN_DGTZ_MallocReadoutBuffer requires char** (ptr must be initialized to NULL)
            ptr = nullptr;
            cmd_local(CAEN_DGTZ_MallocReadoutBuffer(handle, &ptr, &size), "MallocReadoutBuffer");
            if (!ptr || size == 0) {
                // allocation failed or returned empty buffer
                ptr = nullptr;
                size = 0;
                throw std::runtime_error("CAEN: failed to allocate readout buffer");
            }
        }

        ~Buffer() {
            if (ptr) {
                // CAEN_DGTZ_FreeReadoutBuffer takes a char** as per header
                CAEN_DGTZ_FreeReadoutBuffer(&ptr);
                ptr = nullptr;
            }
        }

        Buffer(const Buffer&) = delete;
        Buffer(Buffer&& o) noexcept : ptr(o.ptr), size(o.size), handle(o.handle) {
            o.ptr = nullptr; o.size = 0; o.handle = -1;
        }
        Buffer& operator=(Buffer&& o) noexcept {
            if (this != &o) {
                if (ptr) CAEN_DGTZ_FreeReadoutBuffer(&ptr);
                ptr = o.ptr; size = o.size; handle = o.handle;
                o.ptr = nullptr; o.size = 0; o.handle = -1;
            }
            return *this;
        }

        char* data() const { return ptr; }
        uint32_t bytes() const { return size; }
        [[nodiscard]] bool valid() const { return ptr != nullptr && size > 0; }

    private:
        // small helper to throw with context inside Buffer ctor (can't call Digitizer::cmd)
        void cmd_local(CAEN_DGTZ_ErrorCode ret, const char* where) {
            if (ret != CAEN_DGTZ_Success) {
                char msg[256];
                snprintf(msg, sizeof(msg), "CAEN error %d at %s", static_cast<int>(ret), where);
                throw std::runtime_error(msg);
            }
        }
    };

    Buffer mallocBuffer() const { return Buffer(handle_); }

    // ------------------------------------------------------------------
    // Accessors
    // ------------------------------------------------------------------
    int handle() const noexcept { return handle_; }
    bool isOpen() const noexcept { return handle_ >= 0; }

    CAEN_DGTZ_BoardInfo_t info() const {
        if (!isOpen()) throw std::runtime_error("Digitizer not open");
        CAEN_DGTZ_BoardInfo_t i{};
        cmd(CAEN_DGTZ_GetInfo(handle_, &i), "GetInfo");
        return i;
    }

private:
    void close() noexcept {
        if (handle_ >= 0) {
            // Best-effort close; ignore errors in destructor context
            CAEN_DGTZ_CloseDigitizer(handle_);
            handle_ = -1;
        }
    }

public:

    // Helper struct to hold discovery info (doesn't touch CAEN's struct)
    struct BoardDiscoveryInfo_t {
        CAEN_DGTZ_BoardInfo_t board;  // Original CAEN info
        CAEN_DGTZ_ConnectionType type;  // USB, OpticalLink, etc.
        int linkNum;
        int conetNode;
    };

    // Enumerate all available boards (returns vector of discovery info)
    static std::vector<BoardDiscoveryInfo_t> listConnected() {
        std::vector<BoardDiscoveryInfo_t> list;

        // Try USB devices (linkNum enumerates USB device number; ConetNode param receives the device number)
        for (int usbIdx = 0; usbIdx < 32; ++usbIdx) {
            int handle = -1;
            CAEN_DGTZ_ErrorCode ret = CAEN_DGTZ_OpenDigitizer2(CAEN_DGTZ_USB, nullptr, usbIdx, 0, &handle);
            if (ret == CAEN_DGTZ_Success && handle >= 0) {
                BoardDiscoveryInfo_t info;
                CAEN_DGTZ_GetInfo(handle, &info.board);
                CAEN_DGTZ_CloseDigitizer(handle);

                info.type = CAEN_DGTZ_USB;
                info.linkNum = usbIdx;
                info.conetNode = usbIdx;
                list.emplace_back(std::move(info));
            }
        }

        // Try Optical Link (Conet node scanning)
        for (int linkNum = 0; linkNum < 32; ++linkNum) {
            for (int conetNode = 0; conetNode < 8; ++conetNode) {
                int handle = -1;
                CAEN_DGTZ_ErrorCode ret = CAEN_DGTZ_OpenDigitizer2(CAEN_DGTZ_OpticalLink, nullptr, conetNode, 0, &handle);
                if (ret == CAEN_DGTZ_Success && handle >= 0) {
                    BoardDiscoveryInfo_t info;
                    CAEN_DGTZ_GetInfo(handle, &info.board);
                    CAEN_DGTZ_CloseDigitizer(handle);

                    info.type = CAEN_DGTZ_OpticalLink;
                    info.linkNum = linkNum;
                    info.conetNode = conetNode;
                    list.emplace_back(std::move(info));
                }
            }
        }

        return list;
    }
};
